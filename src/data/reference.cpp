#include "data/bData.hpp"
#include "data/vData.hpp"
#include "data/tokens.hpp"
#include "tools/gtf_data.hpp"
#include "data/reference.hpp"
#include "VarQuin/VarQuin.hpp"
#include "MetaQuin/MetaQuin.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

template <typename Key, typename Value> std::set<Key> getKeys(const std::map<Key, Value> &m)
{
    std::set<Key> keys;

    for(auto i: m)
    {
        keys.insert(i.first);
    }
    
    return keys;
}

/*
 * ------------------------- Transcriptome Analysis -------------------------
 */

template <typename Iter> Base countLocus(const Iter &iter)
{
    Base n = 0;
    
    for (const auto &i : iter)
    {
        n += static_cast<Locus>(i).length();
    }
    
    return n;
}

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
            Standard::addGenomic(i.first);
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
    
    for (const auto i : _mixes)
    {
        // Eg: MixA, MixB etc
        const auto mix = i.first;
        
        // For each of the sequin
        for (const auto j : i.second)
        {
            // Only if it's a validated sequin
            if (_data.count(j.second->id))
            {
                _data.at(j.first).mixes[mix] = j.second->abund;
            }
        }
    }

    assert(!_data.empty());
}

void RnaRef::validate()
{
    auto iIDs = std::set<SequinID>();
    
    for (const auto &i : _impl->gData)
    {
        if (isRnaQuin(i.first))
        {
            iIDs = getKeys(_impl->gData.at(i.first).t2d);
            break;
        }
    }
    
    /*
     * Building rules:
     *
     *   1: Only annoation
     *   2: Only mixture
     *   3: Annotation and mixture
     */
    
    if (_rawMIDs.empty())
    {
        merge(iIDs, iIDs);         // Rule 1
    }
    else if (!iIDs.empty())
    {
        merge(_rawMIDs, iIDs);     // Rule 3
    }
    else
    {
        merge(_rawMIDs, _rawMIDs); // Rule 2
    }
    
    /*
     * Always prefer reference annotation be given. However, if this is not provided, we'll need to
     * work out the RNA structure ourself. Coordinates are not required.
     */
    
    if (_impl->gData.empty())
    {
        for (const auto &i : _rawMIDs)
        {
            TransData t;
            
            t.cID = ChrIS;
            t.tID = i;
            t.gID = Isoform2Gene(i);
            
            GeneData g;
            
            g.cID = ChrIS;
            g.gID = t.gID;
            
            const auto mix = findMix(Mix_1, t.tID);
            assert(mix);
            
            g.l = _impl->gData[ChrIS].g2d[t.gID].l;
            t.l = Locus(1, mix->length);

            // Merge the transcripts...
            g.l.merge(Locus(1, mix->length));
            
            assert(g.l.length() > 1);
            
            _impl->gData[ChrIS].g2d[t.gID] = g;
            _impl->gData[ChrIS].t2d[t.tID] = t;
            _impl->gData[ChrIS].t2g[t.tID] = t.gID;
        }
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

void MetaRef::readBed(const Reader &r)
{
    for (const auto &i : (_impl->bData = readRegions(r)))
    {
        if (!isMetaQuin(i.first))
        {
            Standard::addGenomic(i.first);
        }
    }
}

MC2Intervals MetaRef::mInters() const
{
    return _impl->bData.minters();
}

Base MetaRef::nBaseSyn() const { return _impl->bData.countBaseSyn(isMetaQuin); }
Base MetaRef::nBaseGen() const { return _impl->bData.countBaseGen(isMetaQuin); }

Counts MetaRef::nMicroSyn() const { return _impl->bData.nGeneSyn(isMetaQuin); }
Counts MetaRef::nMicroGen() const { return _impl->bData.nGeneGen(isMetaQuin); }

void MetaRef::validate()
{
    auto bed2ID = [](const BedData &data)
    {
        std::set<SequinID> ids;
        
        std::for_each(data.begin(), data.end(), [&](const std::pair<ChrID, BedChrData> & p)
        {
            if (isMetaQuin(p.first))
            {
                ids.insert(p.first);
            }
        });
        
        A_CHECK(!ids.empty(), "No sequin found in the reference");
        
        return ids;
    };

    if (!_impl->bData.length())
    {
        merge(_rawMIDs, _rawMIDs);
    }
    else if (_rawMIDs.empty())
    {
        merge(bed2ID(_impl->bData));
    }
    else
    {
        merge(_rawMIDs, bed2ID(_impl->bData));
        
        /*
         * Build length for each synthetic genome
         */
        
        for (const auto &i : _impl->bData)
        {
            for (const auto &j : i.second.r2d)
            {
                if (_data.count(i.first))
                {
                    _data.at(i.first).l = j.second.l;
                }
            }
        }
    }
}

/*
 * ------------------------- Variant Analysis -------------------------
 */

struct VarRef::VarRefImpl
{
    // Optimization (no need to use the whole data-structure)
    std::set<SequinID> vIDs, bIDs;
    
    VData vData;
    
    // Information about sequin variants
    std::map<VarKey, SeqVariant> sVars;

    // Regular regions
    BedData bData;
    
    // Trimmed regions
    BedData tData;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

Counts VarRef::nType(Mutation m) const
{
    return _impl->vData.count_(m);
}

Counts VarRef::nGroup(SeqVariant::Group grp) const
{
    return countMap(_impl->sVars, [&](VarKey, const SeqVariant &x)
    {
        return x.group == grp ? 1 : 0;
    });
}

std::set<Variant> VarRef::vars() const
{
    return _impl->vData.vars();
}

const SeqVariant & VarRef::findSeqVar(long key) const
{
    return _impl->sVars.at(key);
}

void VarRef::readGBRef(const Reader &r, Base trim)
{
    _impl->bData = readRegions(r, [&](const ParserBed::Data &x, const ParserProgress &)
    {
        _impl->bIDs.insert(x.name);
    });
    
    RegionOptions o;
    o.trim = trim;

    // Read trimmed regions
    _impl->tData = readRegions(Reader(r), [&](const ParserBed::Data &x, const ParserProgress &) {}, o);

    A_CHECK(!_impl->bData.empty(), "No sampling regions for sampling");
    A_CHECK(!_impl->tData.empty(), "No sampling regions for sampling");
}

void VarRef::readVRef(const Reader &r)
{
    auto throwInvalidRef = [&](const std::string &x)
    {
        throw std::runtime_error(r.src() + " doesn't seem to be a valid VCF reference file. Reason: " + x);
    };
    
    _impl->vData = readVFile(r, [&](const ParserVCF::Data &x, const ParserProgress &)
    {
        A_ASSERT(x.key());
        
        typedef SeqVariant::Group Group;
        
        const auto m1 = std::map<std::string, SeqVariant::Group>
        {
            { "NA12878", Group::NA12878 },
            { "high_gc", Group::HighGC },
            { "long_dinuc_rep", Group::LongDinRep },
            { "long_homopol", Group::LongHompo },
            { "low_gc", Group::LowGC },
            { "long_quadnuc_rep", Group::LongQuadRep },
            { "long_trinuc_rep", Group::LongTrinRep },
            { "short_dinuc_rep", Group::ShortDinRep },
            { "short_hompol", Group::ShortHompo },
            { "short_quadnuc_rep", Group::ShortQuadRep },
            { "very_high_gc", Group::VeryHighGC },
            { "very_low_gc", Group::VeryLowGC },
            { "short_trinuc_rep", Group::ShortTrinRep }
        };
        
        const auto m2 = std::map<std::string, Zygosity>
        {
            { "HOM",     Zygosity::Homozygous  },
            { "HET",     Zygosity::Heterzygous },
            { "HOM_CNV", Zygosity::Homozygous  },
        };
        
        if (!x.opts.count("GR") || !m1.count(x.opts.at("GR")))
        {
            if (!x.opts.count("GR") || !boost::starts_with(x.opts.at("GR"), "COSM"))
            {
                throwInvalidRef("The GR field is not found or invalid");
            }
        }
        
        if (!x.opts.count("ZY") || !m2.count(x.opts.at("ZY")))
        {
            throwInvalidRef("The ZY field is not found or invalid");
        }
        else if (!x.opts.count("CP"))
        {
            throwInvalidRef("The CP field is not found or invalid");
        }
        
        SeqVariant s;

        s.zyg  = m2.at(x.opts.at("ZY"));
        s.copy = stoi(x.opts.at("CP"));
        
        const auto &gr = x.opts.at("GR");
        
        if (boost::starts_with(gr, "COSM"))
        {
            s.info  = gr;
            s.group = Group::Cosmic;
        }
        else
        {
            s.group = m1.at(x.opts.at("GR"));
        }
        
        _impl->vIDs.insert(x.name);

        Concent af;
        
        switch (s.zyg)
        {
            case Zygosity::Somatic:     { af = stoi(x.opts.at("AF")); break; }
            case Zygosity::Homozygous:  { af = 1,0; break; }
            case Zygosity::Heterzygous: { af = 0.5; break; }
        }

        _impl->sVars[x.key()] = s;
        
        add(x.name + "_V", _impl->bData.lengthReg(x.name), af, Mix_1);
        add(x.name + "_R", _impl->bData.lengthReg(x.name), 1-af, Mix_1);
    });
}

C2Intervals VarRef::regions(bool trimmed) const
{
    return trimmed ? _impl->tData.inters() : _impl->bData.inters();
}

MC2Intervals VarRef::mInters() const { return _impl->bData.minters(); }

MergedIntervals<> VarRef::mInters(const ChrID &cID) const
{
    return _impl->bData.minters(cID);
}

bool VarRef::isGermline() const
{
    std::set<Proportion> freqs;
    
    for (const auto &i : _mixes.at(Mix_1))
    {
        freqs.insert(i.second->abund);
    }
    
    return freqs.size() == 3 && freqs.count(0) && freqs.count(0.5) && freqs.count(1);
}

Concent VarRef::findRCon(const SequinID &id) const
{
    return _mixes.at(Mix_1).at(id + "_R")->abund;
}

Concent VarRef::findVCon(const SequinID &id) const
{
    return _mixes.at(Mix_1).at(id + "_V")->abund;
}

Proportion VarRef::findAFreq(const SequinID &id) const
{
    const auto v = findRCon(id);
    const auto r = findRCon(id);
    
    return v / (r + v);
}

bool VarRef::hasAFreq(const SequinID &id) const
{
    return _mixes.at(Mix_1).count(id + "_R");
}

Counts VarRef::countInd(const ChrID &cID) const
{
    return _impl->vData.countInd(cID);
}

Counts VarRef::countIndSyn() const
{
    return _impl->vData.countIndSyn();
}

Counts VarRef::countSNP(const ChrID &cID) const
{
    return _impl->vData.countSNP(cID);
}

Counts VarRef::countSNPSyn() const
{
    return _impl->vData.countSNPSyn();
}

Counts VarRef::nRegs() const { return _impl->bData.count();  }
Counts VarRef::lRegs() const { return _impl->bData.length(); }

void VarRef::validate()
{
    const auto hasVCF = !_impl->vIDs.empty();
    const auto hasBED = !_impl->bIDs.empty();
    
    if (hasVCF && hasBED)
    {
        merge(_impl->vIDs);
    }
    else if (hasBED)
    {
        merge(_impl->bIDs);
    }
    else
    {
        // Try mixture if we don't have anything else
        merge(_rawMIDs, _rawMIDs);
    }
}

const Variant * VarRef::findVar(const ChrID &cID, const Locus &l) const
{
    return _impl->vData.findVar(cID, l);
}

const Variant * VarRef::findVar(const ChrID &cID, long key) const
{
    return _impl->vData.findVar(cID, key);
}
