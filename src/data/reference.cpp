#include "data/tokens.hpp"
#include "tools/bed_data.hpp"
#include "tools/gtf_data.hpp"
#include "tools/vcf_data.hpp"
#include "data/reference.hpp"
#include "VarQuin/VarQuin.hpp"
#include "RnaQuin/RnaQuin.hpp"
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

struct TransRef::TransRefImpl
{
    // Includes synthetic and genome
    GTFData gData;
    
    // Intervals for the genes
    std::map<ChrID, Intervals<>> gInters;
};

TransRef::TransRef() : _impl(new TransRefImpl()) {}

Counts TransRef::countLenSyn() const
{
    return _impl->gData.countLenSyn();
}

Counts TransRef::countLenGen() const
{
    return _impl->gData.countLenGen();
}

Counts TransRef::countUExon(const ChrID &cID) const
{
    return _impl->gData.countUExon(cID);
}

Counts TransRef::countUExonSyn() const
{
    return _impl->gData.countUExonSyn();
}

Counts TransRef::countUExonGen() const
{
    return _impl->gData.countUExonGen();
}

Counts TransRef::countUIntr(const ChrID &cID) const
{
    return _impl->gData.countUIntr(cID);
}

Counts TransRef::countUIntrSyn() const
{
    return _impl->gData.countUIntrSyn();
}

Counts TransRef::countUIntrGen() const
{
    return _impl->gData.countUIntrGen();
}

Counts TransRef::countGeneSyn() const
{
    return _impl->gData.countGeneSyn();
}

Counts TransRef::countGeneGen() const
{
    return _impl->gData.countGeneGen();
}

Counts TransRef::countTransSyn() const
{
    return _impl->gData.countTransSyn();
}

Counts TransRef::countTransGen() const
{
    return _impl->gData.countTransGen();
}

void TransRef::readRef(const Reader &r)
{
    for (const auto &i : (_impl->gData = gtfData(r)))
    {
        if (!Standard::isSynthetic(i.first))
        {
            Standard::addGenomic(i.first);
        }
    }
}

std::map<ChrID, Hist> TransRef::histGene() const
{
    return _impl->gData.histGene();
}

std::map<ChrID, Hist> TransRef::histIsof() const
{
    return _impl->gData.histIsof();
}

Counts TransRef::countGeneSeqs() const
{
    std::set<GeneID> gIDs;
    
    for (const auto &i : _data)
    {
        gIDs.insert(RnaQuin::t2g(i.first));
    }
    
    return gIDs.size();
}

LogFold TransRef::logFoldGene(const GeneID &gID) const
{
    const auto e1 = concent(gID, Mix_1);
    const auto e2 = concent(gID, Mix_2);

    return log2(e2 / e1);
}

LogFold TransRef::logFoldSeq(const IsoformID &iID) const
{
    const auto m = match(iID);
    
    // It's pre-condition that the sequin exists
    assert(m);
    
    const auto e1 = m->concent(Mix_1);
    const auto e2 = m->concent(Mix_2);
    
    return log2(e2 / e1);
}

Concent TransRef::concent(const GeneID &gID, Mixture mix) const
{
    for (const auto &i : _impl->gData)
    {
        if (Standard::isSynthetic(i.first))
        {
            A_ASSERT(!i.second.t2g.empty(), "No transcript found in gene [" + gID + "]");

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

GeneID TransRef::s2g(const SequinID &sID) const
{
    return _impl->gData.at(ChrIS).t2g.at(sID);
}

const TransData *TransRef::findTrans(const ChrID &cID, const TransID &tID) const
{
    if (!_impl->gData.count(cID))
    {
        return nullptr;
    }
    
    assert(!_impl->gData.at(cID).t2d.empty());
    return _impl->gData.at(cID).t2d.count(tID) ? &(_impl->gData.at(cID).t2d[tID]) : nullptr;
}

const GeneData * TransRef::findGene(const ChrID &cID, const GeneID &gID) const
{
    if (!_impl->gData.count(cID))
    {
        return nullptr;
    }
    
    assert(!_impl->gData.at(cID).g2d.empty());
    return _impl->gData.at(cID).g2d.count(gID) ? &(_impl->gData.at(cID).g2d[gID]) : nullptr;
}

std::set<GeneID> TransRef::getGenes(const ChrID &cID) const
{
    std::set<GeneID> ids;
    
    for (const auto &i : _impl->gData.at(cID).g2d)
    {
        ids.insert(i.first);
    }
    
    return ids;
}

std::set<TransID> TransRef::getTrans(const ChrID &cID) const
{
    std::set<GeneID> ids;
    
    for (const auto &i : _impl->gData.at(cID).t2d)
    {
        ids.insert(i.first);
    }
    
    return ids;
}

MergedIntervals<> TransRef::mergedExons(const ChrID &cID) const
{
    return _impl->gData.mergedExons(cID);
}

MC2Intervals TransRef::mergedExons() const
{
    return _impl->gData.mergedExons();
}

MC2Intervals TransRef::meInters(Strand str) const
{
    return _impl->gData.meInters(str);
}

MC2Intervals TransRef::ueInters() const
{
    return _impl->gData.ueInters();
}

MC2Intervals TransRef::uiInters() const
{
    return _impl->gData.uiInters();
}

void TransRef::merge(const std::set<SequinID> &mIDs, const std::set<SequinID> &aIDs)
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
        
        // For each of the mixture defined
        for (const auto j : i.second)
        {
            // Only if it's a validated sequin
            if (_data.count(j.id))
            {
                _data.at(j.id).mixes[mix] = j.abund;
            }
        }
    }

    assert(!_data.empty());
}

void TransRef::validate()
{
    auto iIDs = std::set<SequinID>();
    
    for (const auto &i : _impl->gData)
    {
        if (Standard::isSynthetic(i.first))
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
            t.gID = RnaQuin::t2g(i);
            
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
 * ------------------------- Variant Analysis -------------------------
 */

struct VarRef::VariantPair
{
    const MixtureData *r, *v;
};

struct VarRef::VarRefImpl
{
    // Mixture data
    std::map<Mixture, std::map<SequinID, VariantPair>> data;

    VCFData vData;
    BedData bData;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

void VarRef::readBRef(const Reader &r)
{
    for (const auto &i : (_impl->bData = bedData(r)))
    {
        if (!Standard::isSynthetic(i.first))
        {
            Standard::addGenomic(i.first);
        }
    }
}

void VarRef::readVRef(const Reader &r)
{
    for (const auto &i : (_impl->vData = vcfData(r)))
    {
        if (!Standard::isSynthetic(i.first))
        {
            Standard::addGenomic(i.first);
        }
    }
}

std::map<ChrID, Hist> VarRef::hist() const
{
    return _impl->bData.hist();
}

C2Intervals  VarRef::dInters()    const { return _impl->bData.inters();    }
ID2Intervals VarRef::dIntersSyn() const { return _impl->bData.intersSyn(); }

MC2Intervals VarRef::mInters()  const { return _impl->bData.minters();    }
MC2Intervals VarRef::msInters() const { return _impl->bData.mintersSyn(); }
MC2Intervals VarRef::mgInters() const { return _impl->bData.mintersGen(); }

MergedIntervals<> VarRef::mInters(const ChrID &cID) const
{
    return _impl->bData.minters(cID);
}

bool VarRef::isGermline() const
{
    std::set<Proportion> freqs;
    
    for (const auto &i : _impl->data.at(Mix_1))
    {
        freqs.insert(findAFreq(i.first));
    }
    
    return freqs.size() == 1;
}

Concent VarRef::findRCon(const SequinID &id) const
{
    const auto &p = _impl->data.at(Mix_1).at(baseID(id));
    return p.r->abund;
}

Concent VarRef::findVCon(const SequinID &id) const
{
    const auto &p = _impl->data.at(Mix_1).at(baseID(id));
    return p.v->abund;
}

Proportion VarRef::findAFreq(const SequinID &id) const
{
    const auto &p = _impl->data.at(Mix_1).at(baseID(id));
    const auto &r = p.r;
    const auto &v = p.v;

    return v->abund / (r->abund + v->abund);
}

Fold VarRef::findAFold(const SequinID &id) const
{
    const auto &p = _impl->data.at(Mix_1).at(baseID(id));
    return round(p.r->abund / p.v->abund);
}

Counts VarRef::countInd(const ChrID &cID) const
{
    return _impl->vData.countInd(cID);
}

Counts VarRef::countIndSyn() const
{
    return _impl->vData.countIndSyn();
}

Counts VarRef::countIndGen() const
{
    return _impl->vData.countIndGen();
}
    
Counts VarRef::countSNP(const ChrID &cID) const
{
    return _impl->vData.countSNP(cID);
}

Counts VarRef::countSNPSyn() const
{
    return _impl->vData.countSNPSyn();
}

Counts VarRef::countSNPGen() const
{
    return _impl->vData.countSNPGen();
}

Counts VarRef::countVar() const
{
    return _impl->vData.countVar();
}

Base VarRef::countBaseSyn() const
{
    return _impl->bData.countBaseSyn();
}

Base VarRef::countBaseGen() const
{
    return _impl->bData.countBaseGen();
}

Counts VarRef::countGeneSyn() const
{
    return _impl->bData.countGeneSyn();
}

Counts VarRef::countGeneGen() const
{
    return _impl->bData.countGeneGen();
}

void VarRef::validate()
{
    // Coordinate annotation?
    const auto shouldCoord = !_impl->bData["chr1"].g2d.empty(); // TODO: Need to loop through all chromosomes
    
    if (shouldCoord)
    {
        std::set<SequinID> ids;
        
        for (const auto &i : _impl->bData["chr1"].g2d)
        {
            ids.insert(i.first);
        }
        
        merge(ids);
    }
    else
    {
        throw std::runtime_error("Failed to find VarQuin sequins in the annotation file");
    }
    
    /*
     * Constructing allele frequency for the variants
     */

    for (const auto &i : _mixes)
    {
        const auto &data = _mixes.at(i.first);

        for (const auto &j : data)
        {
            if (isRefID(j.id))
            {
                // Eg: D_1_3_R
                const auto rID = j.id;
                
                // Eg: D_1_3_V
                const auto vID = rID.substr(0, rID.size() - 2) + "_V";
             
                auto rIter = std::find_if(data.begin(), data.end(), [&](const MixtureData &m)
                {
                    return m.id == rID;
                });
                
                auto vIter = std::find_if(data.begin(), data.end(), [&](const MixtureData &m)
                {
                    return m.id == vID;
                });
                
                assert(rIter != data.end() && vIter != data.end());
                
                _impl->data[i.first][baseID(j.id)].r = &(*rIter);
                _impl->data[i.first][baseID(j.id)].v = &(*vIter);
            }
        }
    }
}

const Variant * VarRef::findVar(const ChrID &cID, const Locus &l) const
{
    return _impl->vData.findVar(cID, l);
}

std::map<ChrID, std::map<long, Counts>> VarRef::vHist() const
{
    return _impl->vData.hist();
}

const Variant * VarRef::findVar(const ChrID &cID, long key) const
{
    return _impl->vData.findVar(cID, key);
}