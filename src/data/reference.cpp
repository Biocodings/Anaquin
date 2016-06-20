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
 * ------------------------- Metagenomics Analysis -------------------------
 */

struct MetaRef::MetaRefImpl
{
    std::map<SequinID, Locus> s2l;
};

MetaRef::MetaRef() : _impl(new MetaRefImpl()) {}

void MetaRef::addStand(const SequinID &id, const Locus &l)
{
    _impl->s2l[id] = l;
}

void MetaRef::validate()
{
    /*
     * Validation rule:
     *
     *   1: Standards & Mixtures
     *   2: Mixtures
     *   3: Standards
     */
    
    if (!_rawMIDs.empty() && !_impl->s2l.empty()) // Case 1
    {
        merge(_rawMIDs, getKeys(_impl->s2l));
    }
    else if (!_rawMIDs.empty())                         // Case 2
    {
        merge(_rawMIDs);
    }
    else if (!_impl->s2l.empty())
    {
        merge(getKeys(_impl->s2l));
    }
    else
    {
        throw std::runtime_error("Unknown validation error");
    }

    /*
     * Build length for each synthetic genome
     */
    
    for (const auto &i : _impl->s2l)
    {
        if (_data.count(i.first))
        {
            _data.at(i.first).l = i.second;
        }
    }
}

const SequinData * MetaRef::contains(const GenomeID &id, const Locus &l) const
{
    return _data.count(id) && _data.at(id).l.contains(l) ? &(_data.at(id)) : nullptr;
}

/*
 * ------------------------- Ladder Analysis -------------------------
 */

struct LadderRef::LadderRefImpl
{
    struct JoinData
    {
        const SequinData *A, *B, *C, *D;
        
        inline Concent concent(Mixture m) const
        {
            return A->concent(m) + B->concent(m) + C->concent(m) + D->concent(m);
        }
    };

    std::map<JoinID, JoinData> joined;
};

LadderRef::LadderRef() : _impl(new LadderRefImpl()) {}

Limit LadderRef::limitJoin(const JoinHist &h) const
{
    return Reference<SequinData, DefaultStats>::detectLimit(h, [&](const JoinID &id)
    {
        return &(_impl->joined.at(id));
    });
}

LadderRef::JoinIDs LadderRef::joinIDs() const
{
    JoinIDs ids;

    for (const auto &i : _impl->joined)
    {
        ids.insert(i.first);
    }
    
    return ids;
}

LadderRef::JoinHist LadderRef::joinHist() const
{
    JoinHist h;
    
    for (const auto &i : _impl->joined)
    {
        h[i.first] = 0;
    }

    return h;
}

void LadderRef::validate()
{
    merge(_rawMIDs);
    
    std::vector<std::string> toks;
    
    for (const auto &i : _data)
    {
        // Eg: C_01_A
        Tokens::split(i.first, "_", toks);

        // Eg: C_01
        const auto joinID = toks[0] + "_" + toks[1];
        
        // Eg: A
        const auto typeID = toks[2];
     
        if      (typeID == "A") { _impl->joined[joinID].A = &i.second; }
        else if (typeID == "B") { _impl->joined[joinID].B = &i.second; }
        else if (typeID == "C") { _impl->joined[joinID].C = &i.second; }
        else if (typeID == "D") { _impl->joined[joinID].D = &i.second; }
        else
        {
            throw std::runtime_error("Unknown " + i.first + ". Mixture file of unjoined standards is required.");
        }
    }

    assert(!_impl->joined.empty());
}

void LadderRef::concent(const LadderRef::JoinID &id, Concent &a, Concent &b, Concent &c, Concent &d, Mixture m) const
{
    a = _impl->joined.at(id).A->concent(m);
    b = _impl->joined.at(id).A->concent(m);
    c = _impl->joined.at(id).A->concent(m);
    d = _impl->joined.at(id).A->concent(m);
}

/*
 * ------------------------- Structual Analysis -------------------------
 */

struct StructRef::StructRefImpl
{
    struct Data
    {
        // TODO: Implement me
    };

    std::map<SequinID, Data> data;
};

void StructRef::addStruct(const SequinID &id) const
{
    _impl->data[id];
}

void StructRef::validate()
{
    
}

/*
 * ------------------------- Fusion Analysis -------------------------
 */

struct FusionRef::FusionRefImpl
{
    // Reference junctions in the parent genes
    std::map<SequinID, Locus> juncts;

    // Known fusions
    std::set<KnownFusion> knowns;

    // Standards (fusion genes will have two loci)
    std::map<SequinID, std::vector<Locus>> stands;
    
    // Normal genes in the standards
    std::map<SequinID, Locus> normals;

    // Fusion genes in the standards
    std::map<SequinID, std::vector<Locus>> fusions;

    /*
     * Mapping table, fusion to normal and vice versa
     */

    std::map<SequinID, SequinID> normToFus;
    std::map<SequinID, SequinID> fusToNorm;

    std::map<SequinID, NormalFusion> normFus;
};

FusionRef::FusionRef() : _impl(new FusionRefImpl()) {}

Counts FusionRef::countFusion() const { return _impl->knowns.size(); }
Counts FusionRef::countJuncts() const { return _impl->juncts.size(); }

const FusionRef::NormalFusion * FusionRef::findNormFus(const SequinID &id) const
{
    return _impl->normFus.count(id) ? &_impl->normFus.at(id) : nullptr;
}

void FusionRef::addFusion(const KnownFusion &f)
{
    _impl->knowns.insert(f);
}

void FusionRef::addJunct(const SequinID &id, const Locus &l)
{
    _impl->juncts[id] = l;
}

void FusionRef::addStand(const SequinID &id, const Locus &l)
{
    _impl->stands[id].push_back(l);
}

SequinHist FusionRef::normalHist() const
{
    return createHist(_impl->juncts);
}

SequinHist FusionRef::fusionHist() const
{
    return createHist(_impl->knowns);
}

const FusionRef::KnownFusion * FusionRef::findFusion(const SequinID &id) const
{
    assert(!_impl->knowns.empty());
    
    for (auto &i : _impl->knowns)
    {
        if (i.id == id)
        {
            return &(i);
        }
    }

    return nullptr;
}

const SequinData *FusionRef::findJunct(const Locus &l) const
{
    assert(!_impl->juncts.empty());

    for (auto &i : _impl->juncts)
    {
        if (i.second == l)
        {
            return &(_data.at(i.first));
        }
    }

    return nullptr;
}

void FusionRef::validate()
{
    /*
     * Validation:
     *
     *   1: Mixtures & Fusions (eg: FusExpress)
     *   2: Fusions (eg: FusDiscover)
     *   3: Standards & Mixtures (eg: FusionAlign)
     *   4: Splicing & Fusions & Mixtures (eg: FusDiff)
     *   5: Splicing & Mixtures (eg: FusExpress)
     */

    // Case 4
    if (!_rawMIDs.empty() && !_impl->knowns.empty() && !_impl->juncts.empty())
    {
        if (_impl->knowns.size() != _impl->juncts.size())
        {
            throw std::runtime_error("Number of fusions not equal to splicing. Please check and try again.");
        }

        merge(_rawMIDs);
        
        /*
         * Constructing a mapping between splicing and fusion. While we can't assume the orders, we can
         * assume like: "FG1_12_P2" and "NG1_12_P2".
         */
        
        auto f = [&](const SequinID &x, const SequinID &y)
        {
            return x.substr(2, x.size()-1) == y.substr(2, y.size()-1);
        };
        
        for (auto &i: _impl->knowns)
        {
            for (auto &j : _impl->juncts)
            {
                if (f(i.id, j.first))
                {
                    // Fusion to normal
                    _impl->fusToNorm[i.id] = j.first;
                    
                    // Normal to fusion
                    _impl->normToFus[j.first] = i.id;

                    break;
                }
            }
        }
        
        if (_impl->fusToNorm.size() != _impl->knowns.size())
        {
            throw std::runtime_error("Failed to construct a mapping table. Please check and try again.");
        }
        
        /*
         * Building differential between normal and fusion genes
         */
        
        for (const auto &i : _impl->normToFus)
        {
            NormalFusion d;
            
            const auto normal = match(i.first);
            const auto fusion = match(i.second);

            d.normal = normal->concent(Mix_1);
            d.fusion = fusion->concent(Mix_1);
            
            _impl->normFus[i.first]  = d;
            _impl->normFus[i.second] = d;
        }
    }

    // Case 5
    else if (!_rawMIDs.empty() && !_impl->juncts.empty())
    {
        merge(_rawMIDs, getKeys(_impl->juncts));
    }

    // Case 3
    else if (!_rawMIDs.empty() && !_impl->stands.empty())
    {
        merge(_rawMIDs);

        for (const auto &i : _impl->stands)
        {
            if (_data.count(i.first))
            {
                // Is this a fusion gene?
                if (i.second.size() == 2)
                {
                    _impl->fusions[i.first] = i.second;
                }
                
                // This must be a normal gene
                else
                {
                    _impl->normals[i.first] = i.second[0];
                }
            }
        }

        // That's because fusion genes are formed by breaking of normal genes
        assert(_impl->fusions.size() == _impl->normals.size());
    }
    
    // Case 1
    else if (!_rawMIDs.empty())
    {
        merge(_rawMIDs);
    }
    
    // Case 2
    else if (!_impl->knowns.empty())
    {
        merge(_impl->knowns);
    }

    else
    {
        throw std::runtime_error("Unknown validation");
    }
}

inline bool compare(Base x, Base y, Base fuzzy = 0.0)
{
    return std::abs(x - y) <= fuzzy;
}

const FusionRef::KnownFusion * FusionRef::findFusion(Base x, Base y, Strand o1, Strand o2, double fuzzy) const
{
    for (const auto &f : _impl->knowns)
    {
        // Match in bases?
        const auto b_match = compare(x, f.l1, fuzzy) && compare(y, f.l2, fuzzy);
        
        // Match in orientation?
        const auto s_match = compare(f.s1, f.s1, fuzzy) && compare(f.s2, f.s2, fuzzy);
        
        if (b_match && s_match)
        {
            return &f;
        }
    }
  
    return nullptr;
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

Counts TransRef::countExonSyn() const
{
    return _impl->gData.countExonSyn();
}

Counts TransRef::countExonGen() const
{
    return _impl->gData.countExonGen();
}

Counts TransRef::countIntrSyn() const
{
    return _impl->gData.countIntrSyn();
}

Counts TransRef::countIntrGen() const
{
    return _impl->gData.countIntrGen();
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

Limit TransRef::geneLimit(const SequinHist &hist) const
{
    throw "Not Implemented";
/*
    return Reference<TransData, DefaultStats>::detectLimit(hist, [&](const GeneID &id)
    {
        return findGene(ChrT, id);
    });
*/
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

    _impl->gInters = _impl->gData.gIntervals();
}

std::map<ChrID, Hist> TransRef::histGene() const
{
    return _impl->gData.histGene();
}

Concent TransRef::concent(const GeneID &gID, Mixture m) const
{
    for (const auto &i : _impl->gData)
    {
        if (Standard::isSynthetic(i.first))
        {
            Concent r = 0;
            
            assert(!i.second.t2g.empty());
            
            for (const auto &j : i.second.t2g)
            {
                if (j.second == gID)
                {
                    r += match(j.first)->mixes.at(m);
                }
            }
            
            if (!r)
            {
                throw "Failed to find gene [" + gID + "] for mixture";
            }
            
            return r;
        }
    }
    
    throw "Failed to find gene [" + gID + "] for mixture";
}

const GeneData * TransRef::findGene(const ChrID &cID, const GeneID &gID) const
{
    assert(!_impl->gData.at(cID).g2d.empty());
    return _impl->gData.at(cID).g2d.count(gID) ? &(_impl->gData.at(cID).g2d[gID]) : nullptr;
}

template <typename Iter> const typename Iter::value_type *findList(const Iter &x, const Locus &l, MatchRule m)
{
    for (const auto &i : x)
    {
        if ((m == Exact && i.l == l) || (m == Contains && i.l.contains(l)))
        {
            return &i;
        }
    }
    
    return nullptr;
}

const Interval * TransRef::findGene(const ChrID &cID, const Locus &l, MatchRule m) const
{
    assert(!_impl->gInters.empty());

    switch (m)
    {
        case MatchRule::Exact:    { return _impl->gInters.at(cID).exact(l);    }
        case MatchRule::Overlap:  { return _impl->gInters.at(cID).overlap(l);  }
        case MatchRule::Contains: { return _impl->gInters.at(cID).contains(l); }
    }
}

Intervals<> TransRef::exonInters(const ChrID &cID) const
{
    return _impl->gData.eIntervals(cID);
}

Intervals<> TransRef::intronInters(const ChrID &cID) const
{
    return _impl->gData.iIntervals(cID);
}

SequinHist TransRef::geneHist(const ChrID &cID) const
{
    assert(_impl->gData.count(cID));
    assert(!_impl->gData[cID].g2d.empty());

    return createHist(_impl->gData[cID].g2d);
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
        auto data = TransData();
        
        data.id  = id;
        //data.gID = !_impl->cRaw.rawMapper.empty() ? _impl->cRaw.rawMapper.at(data.id) : "";

        // Add a new entry for the validated sequin
        _data[id] = data;
        
        //assert(!d.id.empty() && !d.gID.empty());
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
    const auto iIDs = _impl->gData.count(ChrT) ? getKeys(_impl->gData.at(ChrT).t2d) : std::set<SequinID>();
    
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
            TransData_ t;
            
            t.cID = ChrT;
            t.tID = i;
            t.gID = RnaQuin::t2g(i);

            GeneData g;
            
            g.cID = ChrT;
            g.gID = t.gID;
            
            _impl->gData[ChrT].g2d[t.gID] = g;
            _impl->gData[ChrT].t2d[t.tID] = t;            
            _impl->gData[ChrT].t2g[t.tID] = t.gID;
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

    /*
     * Data structure for synthetic
     */
    
    // Mixture for bases
    std::map<SequinID, Base> baseMix;

    VCFData vData;
    BedData bData;
    
    std::map<ChrID, Intervals<>> gInters;
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

    _impl->gInters = _impl->bData.gInters();
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

std::map<ChrID, Intervals<>> VarRef::inters() const
{
    return _impl->bData.gInters();
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
    const auto r  = round(p.r->abund / p.v->abund);

    if (r == 2040) { return 2048; }
    if (r == 1019) { return 1024; }
    
    return round(p.r->abund / p.v->abund);
}

Counts VarRef::countInters() const
{
    throw "Not Implemented";
}

Counts VarRef::countSeqs() const
{
    return data().size();
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
    /*
     * Rules:
     *
     *   1: Annotation (eg: VarAlign and VarDiscover)
     *   2: Mixture (eg: VarAlign)
     *   3: Variants & mixture (eg: VarFrequency)
     */
    
    // Rule: 2 and 3
    if (!_rawMIDs.empty())
    {
        merge(_rawMIDs);
    }
    
    // Rule: 1
    else if (!_impl->bData[ChrT].g2d.empty())
    {
        std::set<SequinID> ids;
        
        for (const auto &i : _impl->bData[ChrT].g2d)
        {
            ids.insert(i.first);
        }

        merge(ids);
    }

    else
    {
        throw std::runtime_error("Failed to validate for VarQuin");
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
    
    /*
     * Building the input concentration for the sequins
     */
    
    if (_mixes.empty())
    {
        for (const auto &i : _data)
        {
            _impl->baseMix[baseID(i.first)].id = baseID(i.first);
        }
    }
    else
    {
        for (const auto &i : _mixes)
        {
            for (const auto &j : _mixes.at(i.first))
            {
                _impl->baseMix[baseID(j.id)].id = baseID(j.id);
                _impl->baseMix[baseID(j.id)].total[i.first] += j.abund;
            }
        }
    }
    
    assert(!_impl->baseMix.empty());
}

const VarRef::Base * VarRef::findGene(const SequinID &id, Mixture mix) const
{
    throw "Not Implemented";
    //return _impl->baseMix.count(id) ? &(_impl->baseMix.at(id)) : nullptr;
}

const Variant * VarRef::findVar(const ChrID &cID, const Locus &l) const
{
    return _impl->vData.findVar(cID, l);
}

std::map<ChrID, Intervals<>> VarRef::intersGen() const
{
    return _impl->bData.gIntersGen();
}

std::map<ChrID, std::map<long, Counts>> VarRef::vHist() const
{
    return _impl->vData.hist();
}

Interval * VarRef::findGeno(const ChrID &cID, const Locus &l, MatchRule rule) const
{
    throw "NOT IMPLEMENTED";
    
//    assert(!Standard::isSynthetic(cID));
//    
//    if (!_impl->genome.count(cID))
//    {
//        return nullptr;
//    }
//    
//    switch (rule)
//    {
//        case MatchRule::Contains: { return _impl->genome.at(cID).contains(l); }
//        case MatchRule::Overlap:  { return _impl->genome.at(cID).overlap(l);  }
//        case MatchRule::Exact:
//        {
//            throw "Not Implemented";
//        }
//    }
}

Interval * VarRef::findGeno(const Locus &l) const
{
    throw "NOT IMPLEMENTED";
//    return _impl->genome.at(genoID()).contains(l);
}

const Variant * VarRef::findVar(const ChrID &cID, long key) const
{
    return _impl->vData.findVar(cID, key);
}

Counts VarRef::countIntervals(const ChrID &cID) const
{
    throw "Not Implemented";
//    return _impl->genome.at(cID).size();
}