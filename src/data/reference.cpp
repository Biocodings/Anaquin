#include "data/bData.hpp"
#include "data/vData.hpp"
#include "tools/tools.hpp"
#include "data/tokens.hpp"
#include "tools/gtf_data.hpp"
#include "data/reference.hpp"
#include "parsers/parser_vcf.hpp"

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
    // Empty Implementation
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
    // Empty Implementation
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

Counts VarRef::nSRep() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isShortRepeat() ? 1 : 0;
    });
}

Counts VarRef::nLRep() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isLongRepeat() ? 1 : 0;
    });
}

Counts VarRef::nLowGC() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isLowGC() ? 1 : 0;
    });
}

Counts VarRef::nHighGC() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isHighGC() ? 1 : 0;
    });
}

Counts VarRef::nVLowGC() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isVLowGC() ? 1 : 0;
    });
}

Counts VarRef::nVHighGC() const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.isVHighGC() ? 1 : 0;
    });
}

Base VarRef::len(Genotype x) const
{
    Base n = 0;

    for (const auto &i : _v1->sVars)
    {
        if (i.second.gt == x)
        {
            for (const auto &j : _v1->data.vars())
            {
                if (j.key() == i.first)
                {
                    n += _r2->len(j.cID, j.name);
                    break;
                }
            }
        }
    }

    return n;
}

Base VarRef::len(Variation x) const
{
    Base n = 0;

    for (const auto &i : _v1->data.vars())
    {
        if (i.type() == x)
        {
            n += _r2->len(i.cID, i.name);
        }
    }
    
    return n;
}

Base VarRef::len(SequinVariant::Context x) const
{
    Base n = 0;

    for (const auto &i : _v1->sVars)
    {
        if (i.second.ctx == x)
        {
            for (const auto &j : _v1->data.vars())
            {
                if (j.key() == i.first)
                {
                    n += _r2->len(j.cID, j.name);
                    break;
                }
            }
        }
    }

    return n;
}

Counts VarRef::nGeno1(Genotype g) const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.gt == g ? 1 : 0;
    });
}

Counts VarRef::nType1(Variation x) const
{
    return _v1->data.count_(x);
}

SequinVariant::Context VarRef::ctx1(const Variant &x) const
{
    return _v1->sVars.at(x.key()).ctx;
}

Counts VarRef::nCtx1(SequinVariant::Context c) const
{
    return countMap(_v1->sVars, [&](VarKey, const SequinVariant &x)
    {
        return x.ctx == c ? 1 : 0;
    });
}

std::set<Variant> VarRef::v1() const
{
    return _v1->data.vars();
}

const SequinVariant * VarRef::findVar(const SequinID &x) const
{
    for (const auto &i : _v1->data.vars())
    {
        if (i.name == x)
        {
            return &(findSeqVar1(i.key()));
        }
    }

    return nullptr;
}

const SequinVariant & VarRef::findSeqVar1(long key) const
{
    return _v1->sVars.at(key);
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
        case Tool::VarKmer:
        case Tool::VarKStats:
        {
            build(r.l1);
            break;
        }

        case Tool::VarCopy:
        {
            const auto inter = intersect(r.r1->names(), r.l1->seqs);
            merge(inter.inters);
            
            filter(r.l1, inter.diffs);
            filter(r.r1, inter.diffs);
            filter(r.r2, inter.diffs);

            build(r.l1, r.r1, r.r2);
            break;
        }

        case Tool::VarFlip:
        {
            merge(r.r1->names());
            build(r.r1);
            break;
        }
            
        case Tool::VarTrim:
        case Tool::VarAlign:
        {
            merge(r.r1->names());
            build(r.r1, r.r2);
            break;
        }

        case Tool::VarCalibrate:
        {
            merge(r.r1->names());
            build(r.r1, r.r2);
            break;
        }
            
        case Tool::VarPartition:
        {
            merge(r.r1->names());
            build(r.r1, r.r2, r.r3);
            break;
        }

        case Tool::VarConjoint:
        {
            build(r.l1, r.l2, r.t1, r.t2);
            break;
        }

        case Tool::VarMutation:
        case Tool::VarStructure:
        {
            merge(r.v1->vIDs);
            build(r.r1, r.r2, r.v1, r.v2);
            break;
        }

        default : { break; }
    }
}

const Variant * VarRef::findV1(const ChrID &x, const Locus &l) const
{
    return _v1->data.findVar(x, l);
}

const Variant * VarRef::findV2(const ChrID &x, const Locus &l) const
{
    return _v2->data.findVar(x, l);
}
