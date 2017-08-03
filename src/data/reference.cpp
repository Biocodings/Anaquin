#include "data/bData.hpp"
#include "data/vData.hpp"
#include "tools/tools.hpp"
#include "data/tokens.hpp"
#include "tools/gtf_data.hpp"
#include "data/reference.hpp"
#include "VarQuin/VarQuin.hpp"
#include "MetaQuin/MetaQuin.hpp"
#include "parsers/parser_vcf.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

struct IntersectResults
{
    std::set<SequinID> diffs, inters;
};

template <typename T> IntersectResults intersect(const std::set<T> &t1, const std::set<T> &t2)
{
    std::set<SequinID> x, y;
    
    for (const auto &i : t1) { x.insert(static_cast<SequinID>(i)); }
    for (const auto &i : t2) { y.insert(static_cast<SequinID>(i)); }
    
    A_ASSERT(!x.empty() && !y.empty());
    
    IntersectResults c;
    
    std::set_intersection(x.begin(),
                          x.end(),
                          y.begin(),
                          y.end(),
                          std::inserter(c.inters, c.inters.begin()));

    std::set_difference(x.begin(),
                        x.end(),
                        y.begin(),
                        y.end(),
                        std::inserter(c.diffs, c.diffs.begin()));

    std::set_difference(y.begin(),
                        y.end(),
                        x.begin(),
                        x.end(),
                        std::inserter(c.diffs, c.diffs.begin()));

    return c;
}

/*
 * ------------------------- Transcriptome Analysis -------------------------
 */

struct RnaRef::RnaRefImpl
{
    // Empty Implementation
};

RnaRef::RnaRef() : _impl(new RnaRefImpl()) {}

void RnaRef::validate(Tool x, const UserReference &r)
{
    switch (x)
    {
        case Tool::RnaAlign:      { build(r.g1); break; }
        case Tool::RnaFoldChange:
        case Tool::RnaExpress:    { build(r.l1, r.l2, r.l3, r.l4, r.l5, r.l6);       break; }
        case Tool::RnaAssembly:   { build(r.l1, r.l2, r.l3, r.l4, r.l5, r.l6, r.g1); break; }
        default:                  { break; }
    }
}

/*
 * ------------------------- Metagenomic Analysis -------------------------
 */

struct MetaRef::MetaRefImpl
{
    BedData bData;
};

MetaRef::MetaRef() : _impl(new MetaRefImpl()) {}

void MetaRef::validate(Tool x, const UserReference &r)
{
    switch (x)
    {
        case Tool::MetaCoverage:
        case Tool::MetaAssembly:  { build(r.l1, r.r1); break; }
        case Tool::MetaSubsample: { build(r.r1);       break; }
        default: { break; }
    }    
}

/*
 * ------------------------- Variant Analysis -------------------------
 */

struct VarRef::VarRefImpl
{
    // Optimization (no need to use the whole data-structure)
    std::set<SequinID> vIDs;
    
    VData vData;
    
    // Information about sequin variants
    std::map<VarKey, SequinVariant> sVars;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

Counts VarRef::nCNV(int c) const
{
    return countMap(_impl->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.copy == c ? 1 : 0;
    });
}

Counts VarRef::nGeno(Genotype g) const
{
    return countMap(_impl->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.gt == g ? 1 : 0;
    });
}

Counts VarRef::nType(Variation x) const
{
    return _impl->vData.count_(x);
}

Counts VarRef::nContext(SequinVariant::Context c) const
{
    return countMap(_impl->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.ctx == c ? 1 : 0;
    });
}

std::set<Variant> VarRef::vars() const
{
    return _impl->vData.vars();
}

const SequinVariant & VarRef::findSeqVar(long key) const
{
    return _impl->sVars.at(key);
}

template <typename F, typename T = VData> T parseVCF2(const Reader &r, F f)
{
    T t;
    
    ParserVCF::parse(r, [&](const Variant &x)
    {
        t[x.cID].b2v[x.l.start] = x;
        t[x.cID].m2v[x.type()].insert(x);
        f(x);
    });
    
    return t;
}

void VarRef::readVRef(const Reader &r)
{
    typedef SequinVariant::Context Context;
    
    _impl->vData = parseVCF2(r, [&](const Variant &x)
    {
        A_ASSERT(x.key());

        auto longVar = [&]()
        {
            SequinVariant s;
            
            s.gt = Genotype::Heterzygous;

            _impl->vIDs.insert(x.name);
            _impl->sVars[x.key()] = s;
        };
        
        auto shortVar = [&]()
        {
            const auto m1 = std::map<std::string, SequinVariant::Context>
            {
                { "cancer",     Context::Cancer       },
                { "common",     Context::Common       },
                { "high_gc",    Context::HighGC       },
                { "long_di",    Context::LongDinRep   },
                { "long_homo",  Context::LongHompo    },
                { "low_gc",     Context::LowGC        },
                { "long_quad",  Context::LongQuadRep  },
                { "long_tri",   Context::LongTrinRep  },
                { "short_di",   Context::ShortDinRep  },
                { "short_homo", Context::ShortHompo   },
                { "short_quad", Context::ShortQuadRep },
                { "v_high_gc",  Context::VeryHighGC   },
                { "v_low_gc",   Context::VeryLowGC    },
                { "short_tri",  Context::ShortTrinRep }
            };
            
            const auto m2 = std::map<std::string, Genotype>
            {
                { "SOM",     Genotype::Somatic     },
                { "HOM",     Genotype::Homozygous  },
                { "HOM_CNV", Genotype::Homozygous  },
                { "HET",     Genotype::Heterzygous },
            };
            
            auto throwInvalidRef = [&](const std::string &x)
            {
                throw std::runtime_error(r.src() + " doesn't seem to be a valid VCF reference file. Reason: " + x);
            };

            if (!x.ifs.count("CX") || !m1.count(x.ifs.at("CX")))
            {
                throwInvalidRef("The CX field is not found or invalid");
            }
            else if (!x.ifs.count("GT") || !m2.count(x.ifs.at("GT")))
            {
                throwInvalidRef("The GT field is not found or invalid");
            }
            
            SequinVariant s;
            
            s.gt   = m2.at(x.ifs.at("GT"));
            s.ctx  = m1.at(x.ifs.at("CX"));
            s.copy = x.iff.at("CP");
            
            _impl->vIDs.insert(x.name);
            
            Concent af;
            
            switch (s.gt)
            {
                case Genotype::Somatic:     { af = x.allF; break; }
                case Genotype::Homozygous:  { af = 1,0;    break; }
                case Genotype::Heterzygous: { af = 0.5;    break; }
            }
            
            _impl->sVars[x.key()] = s;
        };
        
        if (x.isSV()) { longVar();  }
        else          { shortVar(); }
    });
}

/*
 * Filter out all reference sequin regions
 */

static void filter(std::shared_ptr<BedData> x, const std::set<SequinID> &ids)
{
    for (const auto &i : ids)
    {
        for (auto &j : *x)
        {
            if (j.second.r2d.count(i))
            {
                j.second.r2d.erase(i);
            }
        }
    }
    
    for (auto i = x->cbegin(); i != x->cend();)
    {
        if (i->second.r2d.empty())
        {
            i = x->erase(i);
        }
        else
        {
            ++i;
        }
    }
}

/*
 * Filter out all reference sequins in ladder
 */

static void filter(std::shared_ptr<Ladder> x, const std::set<SequinID> &ids)
{
    for (const auto &id : ids)
    {
        x->remove(id);
    }
}

void VarRef::validate(Tool x, const UserReference &r)
{
    switch (x)
    {
        case Tool::VarKSomatic:
        {
            build(r.l1, r.l2, r.l3, r.l4, r.l5);
            break;
        }

        case Tool::VarCopy:
        {
            const auto inter = intersect(r.r1->seqs(), r.l1->seqs);

            merge(inter.inters);
            
            filter(r.l1, inter.diffs);
            filter(r.r1, inter.diffs);
            filter(r.r2, inter.diffs);

            build(r.l1, r.r1, r.r2);
            break;
        }

        case Tool::VarTrim:
        case Tool::VarAlign:
        {
            merge(r.r1->seqs());
            build(r.r1, r.r2);
            break;
        }

        case Tool::VarSample:
        {
            merge(r.r1->seqs());
            build(r.r1, r.r2);
            break;
        }

        case Tool::VarConjoint:
        {
            build(r.l1, r.l2, r.t1, r.t2);
            break;
        }

        case Tool::VarDetect:
        case Tool::VarSomatic:
        case Tool::VarStructure:
        {
            merge(_impl->vIDs);
            build(r.r1, r.r2);
            break;
        }

        default : { break; }
    }
}

const Variant * VarRef::findVar(const ChrID &id, const Locus &l) const
{
    return _impl->vData.findVar(id, l);
}
