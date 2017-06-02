#include "data/bData.hpp"
#include "data/vData.hpp"
#include "tools/tools.hpp"
#include "data/tokens.hpp"
#include "tools/gtf_data.hpp"
#include "data/reference.hpp"
#include "VarQuin/VarQuin.hpp"
#include "MetaQuin/MetaQuin.hpp"
#include "parsers/parser_vcf2.hpp"
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
    // Includes synthetic and genome
    GTFData gData;
    
    // Intervals for the genes
    std::map<ChrID, Intervals<>> gInters;
};

RnaRef::RnaRef() : _impl(new RnaRefImpl()) {}

Counts RnaRef::countLenSyn() const
{
    return _impl->gData.countLenSyn();
}

Counts RnaRef::countLenGen() const
{
    return _impl->gData.countLenGen();
}

Counts RnaRef::countUExon(const ChrID &cID) const
{
    return _impl->gData.countUExon(cID);
}

Counts RnaRef::countUExonSyn() const
{
    return _impl->gData.countUExonSyn();
}

Counts RnaRef::countUExonGen() const
{
    return _impl->gData.countUExonGen();
}

Counts RnaRef::countUIntr(const ChrID &cID) const
{
    return _impl->gData.countUIntr(cID);
}

Counts RnaRef::countUIntrSyn() const
{
    return _impl->gData.countUIntrSyn();
}

Counts RnaRef::countUIntrGen() const
{
    return _impl->gData.countUIntrGen();
}

Counts RnaRef::nGeneSyn() const
{
    return _impl->gData.nGeneSyn();
}

Counts RnaRef::nGeneGen() const
{
    return _impl->gData.nGeneGen();
}

Counts RnaRef::countTransSyn() const
{
    return _impl->gData.countTransSyn();
}

Counts RnaRef::countTransGen() const
{
    return _impl->gData.countTransGen();
}

void RnaRef::readRef(const Reader &r)
{
    for (const auto &i : (_impl->gData = gtfData(r)))
    {
        if (!isRnaQuin(i.first))
        {
//            Standard::addGenomic(i.first);
        }
    }
}

std::map<ChrID, Hist> RnaRef::histGene() const
{
    return _impl->gData.histGene();
}

std::map<ChrID, Hist> RnaRef::histIsof() const
{
    return _impl->gData.histIsof();
}

Counts RnaRef::nGeneSeqs() const
{
    std::set<GeneID> gIDs;
    
    for (const auto &i : _data)
    {
        gIDs.insert(Isoform2Gene(i.first));
    }
    
    return gIDs.size();
}

LogFold RnaRef::logFoldGene(const GeneID &gID) const
{
    const auto e1 = concent(gID, Mix_1);
    const auto e2 = concent(gID, Mix_2);

    return log2(e2 / e1);
}

LogFold RnaRef::logFoldSeq(const IsoformID &iID) const
{
    const auto m = match(iID);
    
    // It's pre-condition that the sequin exists
    assert(m);
    
    const auto e1 = m->concent(Mix_1);
    const auto e2 = m->concent(Mix_2);
    
    return log2(e2 / e1);
}

Concent RnaRef::concent(const GeneID &gID, Mixture mix) const
{
    for (const auto &i : _impl->gData)
    {
        if (isRnaQuin(i.first))
        {
            A_CHECK(!i.second.t2g.empty(), "No transcript found in gene [" + gID + "]");

            Concent r = 0;

            for (const auto &j : i.second.t2g)
            {
                // Does this transcript belong to the gene?
                if (j.second == gID)
                {
                    // Add up the concentration
                    r += match(j.first)->mixes.at(mix);
                }
            }

            if (!r)
            {
                A_THROW("Concentration is zero for gene [" + gID + "]");
            }

            return r;
        }
    }
    
    A_THROW("Failed to find gene [" + gID + "]");

    // Never executed
    return Concent();
}

GeneID RnaRef::s2g(const SequinID &sID) const
{
    return _impl->gData.at(ChrIS).t2g.at(sID);
}

const TransData *RnaRef::findTrans(const ChrID &cID, const TransID &tID) const
{
    if (!_impl->gData.count(cID))
    {
        return nullptr;
    }
    
    assert(!_impl->gData.at(cID).t2d.empty());
    return _impl->gData.at(cID).t2d.count(tID) ? &(_impl->gData.at(cID).t2d[tID]) : nullptr;
}

const GeneData * RnaRef::findGene(const ChrID &cID, const GeneID &gID) const
{
    if (!_impl->gData.count(cID))
    {
        return nullptr;
    }
    
    assert(!_impl->gData.at(cID).g2d.empty());
    return _impl->gData.at(cID).g2d.count(gID) ? &(_impl->gData.at(cID).g2d[gID]) : nullptr;
}

std::set<GeneID> RnaRef::getGenes(const ChrID &cID) const
{
    std::set<GeneID> ids;
    
    for (const auto &i : _impl->gData.at(cID).g2d)
    {
        ids.insert(i.first);
    }
    
    return ids;
}

std::set<TransID> RnaRef::getTrans(const ChrID &cID) const
{
    std::set<GeneID> ids;
    
    for (const auto &i : _impl->gData.at(cID).t2d)
    {
        ids.insert(i.first);
    }
    
    return ids;
}

MergedIntervals<> RnaRef::mergedExons(const ChrID &cID) const
{
    return _impl->gData.mergedExons(cID);
}

MC2Intervals RnaRef::mergedExons() const
{
    return _impl->gData.mergedExons();
}

MC2Intervals RnaRef::meInters(Strand str) const
{
    return _impl->gData.meInters(str);
}

MC2Intervals RnaRef::ueInters() const
{
    return _impl->gData.ueInters();
}

MC2Intervals RnaRef::uiInters() const
{
    return _impl->gData.uiInters();
}

void RnaRef::merge(const std::set<SequinID> &mIDs, const std::set<SequinID> &aIDs)
{
    assert(!mIDs.empty() && !aIDs.empty());
    
    std::vector<SequinID> diffs, inters;
    
    /*
     * Check for any sequin defined in mixture but not in annotation
     */
    
    std::set_difference(mIDs.begin(),
                        mIDs.end(),
                        aIDs.begin(),
                        aIDs.end(),
                        std::back_inserter(diffs));
    
    /*
     * Check for any sequin defined in both mixture and annotation
     */
    
    std::set_intersection(mIDs.begin(),
                          mIDs.end(),
                          aIDs.begin(),
                          aIDs.end(),
                          std::back_inserter(inters));
    
    /*
     * Construct a set of validated sequins. A valid sequin is one in which it's
     * defined in both mixture and annoation.
     */
    
    std::for_each(inters.begin(), inters.end(), [&](const SequinID &id)
    {
        auto data = SequinData();
        
        data.id = id;

        // Add a new entry for the validated sequin
        _data[id] = data;

        assert(!id.empty());
    });
    
    /*
     * Now, we have a list of validated sequins. Use those sequins to combine information
     * from mixtures and annotations.
     */
    
//    for (const auto i : _mixes)
//    {
//        // Eg: MixA, MixB etc
//        const auto mix = i.first;
//        
//        // For each of the sequin
//        for (const auto j : i.second)
//        {
//            // Only if it's a validated sequin
//            if (_data.count(j.second->id))
//            {
//                _data.at(j.first).mixes[mix] = j.second->abund;
//            }
//        }
//    }

    assert(!_data.empty());
}

void RnaRef::validate(Tool, const UserReference &r)
{
//    auto iIDs = std::set<SequinID>();
//    
//    for (const auto &i : _impl->gData)
//    {
//        if (isRnaQuin(i.first))
//        {
//            iIDs = keys(_impl->gData.at(i.first).t2d);
//            break;
//        }
//    }
//    
//    /*
//     * Building rules:
//     *
//     *   1: Only annoation
//     *   2: Only mixture
//     *   3: Annotation and mixture
//     */
//    
//    if (_rawMIDs.empty())
//    {
//        merge(iIDs, iIDs);         // Rule 1
//    }
//    else if (!iIDs.empty())
//    {
//        merge(_rawMIDs, iIDs);     // Rule 3
//    }
//    else
//    {
//        merge(_rawMIDs, _rawMIDs); // Rule 2
//    }
//    
//    /*
//     * Always prefer reference annotation be given. However, if this is not provided, we'll need to
//     * work out the RNA structure ourself. Coordinates are not required.
//     */
//    
//    if (_impl->gData.empty())
//    {
//        for (const auto &i : _rawMIDs)
//        {
//            TransData t;
//            
//            t.cID = ChrIS;
//            t.tID = i;
//            t.gID = Isoform2Gene(i);
//            
//            GeneData g;
//            
//            g.cID = ChrIS;
//            g.gID = t.gID;
//            
//            const auto mix = findMix(Mix_1, t.tID);
//            assert(mix);
//            
//            g.l = _impl->gData[ChrIS].g2d[t.gID].l;
//            t.l = Locus(1, mix->length);
//
//            // Merge the transcripts...
//            g.l.merge(Locus(1, mix->length));
//            
//            assert(g.l.length() > 1);
//            
//            _impl->gData[ChrIS].g2d[t.gID] = g;
//            _impl->gData[ChrIS].t2d[t.tID] = t;
//            _impl->gData[ChrIS].t2g[t.tID] = t.gID;
//        }
//    }
}

/*
 * ------------------------- Metagenomic Analysis -------------------------
 */

struct MetaRef::MetaRefImpl
{
    BedData bData;
};

MetaRef::MetaRef() : _impl(new MetaRefImpl()) {}

void MetaRef::readBed(const Reader &r)
{
    for (const auto &i : (_impl->bData = readRegions(r)))
    {
        if (!isMetaQuin(i.first))
        {
//            Standard::addGenomic(i.first);
        }
    }
}

MC2Intervals MetaRef::mInters() const
{
    return _impl->bData.minters();
}

Base MetaRef::nBaseSyn() const { return 0; /*return _impl->bData.countBaseSyn(isMetaQuin);*/ }
Base MetaRef::nBaseGen() const { return 0; /*return _impl->bData.countBaseGen(isMetaQuin);*/ }

Counts MetaRef::nMicroSyn() const { return _impl->bData.nGeneSyn(isMetaQuin); }
Counts MetaRef::nMicroGen() const { return _impl->bData.nGeneGen(isMetaQuin); }

void MetaRef::validate(Tool, const UserReference &r)
{
//    auto bed2ID = [](const BedData &data)
//    {
//        std::set<SequinID> ids;
//        
//        std::for_each(data.begin(), data.end(), [&](const std::pair<ChrID, BedChrData> & p)
//        {
//            if (isMetaQuin(p.first))
//            {
//                ids.insert(p.first);
//            }
//        });
//        
//        A_CHECK(!ids.empty(), "No sequin found in the reference");
//        
//        return ids;
//    };
//
//    if (!_impl->bData.length())
//    {
//        merge(_rawMIDs, _rawMIDs);
//    }
//    else if (_rawMIDs.empty())
//    {
//        merge(bed2ID(_impl->bData));
//    }
//    else
//    {
//        merge(_rawMIDs, bed2ID(_impl->bData));
//        
//        /*
//         * Build length for each synthetic genome
//         */
//        
//        for (const auto &i : _impl->bData)
//        {
//            for (const auto &j : i.second.r2d)
//            {
//                if (_data.count(i.first))
//                {
//                    _data.at(i.first).l = j.second.l;
//                }
//            }
//        }
//    }
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

            if (!x.opts.count("CX") || !m1.count(x.opts.at("CX")))
            {
                throwInvalidRef("The CX field is not found or invalid");
            }
            else if (!x.opts.count("GT") || !m2.count(x.opts.at("GT")))
            {
                throwInvalidRef("The GT field is not found or invalid");
            }
            
            SequinVariant s;
            
            s.gt   = m2.at(x.opts.at("GT"));
            s.ctx  = m1.at(x.opts.at("CX"));
            s.copy = stoi(x.opts.at("CP"));
            
            _impl->vIDs.insert(x.name);
            
            Concent af;
            
            switch (s.gt)
            {
                case Genotype::Somatic:     { af = x.allF; break; }
                case Genotype::Homozygous:  { af = 1,0;    break; }
                case Genotype::Heterzygous: { af = 0.5;    break; }
            }
            
            _impl->sVars[x.key()] = s;
            _mixes[Mix_1][x.name] = std::shared_ptr<MixtureData>(new MixtureData(x.name, 1000, af));
        };
        
        if (x.isSV()) { longVar();  }
        else          { shortVar(); }
    });
}

Proportion VarRef::findAFreq(const SequinID &x) const
{
    return _mixes.at(Mix_1).at(x)->abund;
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
            build(r.r1);
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
            build(r.l1, r.l2);
            break;
        }

        case Tool::VarDetect:
        case Tool::VarStructure:
        {
            merge(_impl->vIDs);
            build(r.r1);
            break;
        }

        default : { break; }
    }
    
    A_ASSERT(!seqs().empty() || !l1Seqs().empty());
}

const Variant * VarRef::findVar(const ChrID &id, const Locus &l) const
{
    return _impl->vData.findVar(id, l);
}
