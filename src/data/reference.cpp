#include <boost/format.hpp>
#include "data/tokens.hpp"
#include "data/reference.hpp"
#include "VarQuin/VarQuin.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

template <typename T> SequinHist createHist(const T& t)
{
    SequinHist hist;
    
    for (const auto &i : t)
    {
        hist[i.first] = 0;
    }
    
    return hist;
}

template <typename T> SequinHist createHist(const std::set<T> & t)
{
    SequinHist hist;
    
    for (const auto &i : t)
    {
        hist[i] = 0;
    }
    
    assert(!hist.empty());

    return hist;
}

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
    /*
     * Raw data
     */

    std::map<SequinID, Locus> rawStands;
};

MetaRef::MetaRef() : _impl(new MetaRefImpl()) {}

void MetaRef::addStand(const SequinID &id, const Locus &l)
{
    _impl->rawStands[id] = l;
}

void MetaRef::validate()
{
    /*
     * Validation rule:
     *
     *   1: Standards & Mixtures (eg: MetaAlign)
     *   2: Mixtures  (eg: MetaAssembly, MetaExpress)
     *   3: Standards (eg: MetaCoverage)
     */
    
    if (!_rawMIDs.empty() && !_impl->rawStands.empty()) // Case 1
    {
        merge(_rawMIDs, getKeys(_impl->rawStands));
    }
    else if (!_rawMIDs.empty())                         // Case 2
    {
        merge(_rawMIDs);
    }
    else if (!_impl->rawStands.empty())
    {
        merge(getKeys(_impl->rawStands));
    }
    else
    {
        throw std::runtime_error("Unknown validation error");
    }

    for (const auto &i : _impl->rawStands)
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
        
        inline Concent abund(Mixture m) const
        {
            return A->abund(m) + B->abund(m) + C->abund(m) + D->abund(m);
        }
    };

    std::map<JoinID, JoinData> joined;
};

LadderRef::LadderRef() : _impl(new LadderRefImpl()) {}

Limit LadderRef::limitJoin(const JoinHist &h) const
{
    return Reference<SequinData, DefaultStats>::absolute(h, [&](const JoinID &id)
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

void LadderRef::abund(const LadderRef::JoinID &id, Concent &a, Concent &b, Concent &c, Concent &d, Mixture m) const
{
    a = _impl->joined.at(id).A->abund(m);
    b = _impl->joined.at(id).A->abund(m);
    c = _impl->joined.at(id).A->abund(m);
    d = _impl->joined.at(id).A->abund(m);
}

/*
 * ------------------------- Fusion Analysis -------------------------
 */

struct FusionRef::FusionRefImpl
{
    // Known splices
    std::map<SequinID, Locus> splices;

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

    std::map<SequinID, SpliceChimeric> spliceChim;
};

FusionRef::FusionRef() : _impl(new FusionRefImpl()) {}

Counts FusionRef::countFusion() const { return _impl->knowns.size();  }
Counts FusionRef::countSplice() const { return _impl->splices.size(); }

const FusionRef::SpliceChimeric * FusionRef::findSpliceChim(const SequinID &id) const
{
    return _impl->spliceChim.count(id) ? &_impl->spliceChim.at(id) : nullptr;
}

SequinID FusionRef::normalToFusion(const SequinID &id) const
{
    return _impl->normToFus.at(id);
}

void FusionRef::addFusion(const KnownFusion &f)
{
    _impl->knowns.insert(f);
}

void FusionRef::addSplice(const SequinID &id, const Locus &l)
{
    _impl->splices[id] = l;
}

void FusionRef::addStand(const SequinID &id, const Locus &l)
{
    _impl->stands[id].push_back(l);
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

const SequinData *FusionRef::findSplice(const Locus &l) const
{
    assert(!_impl->splices.empty());

    for (auto &i : _impl->splices)
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
     *   1: Mixtures & Fusions (eg: FusionExpress)
     *   2: Fusions (eg: FusionDiscover)
     *   3: Standards & Mixtures (eg: FusionAlign)
     *   4: Splicing & Fusions & Mixtures (eg: FusionDiff)
     *   5: Splicing & Mixtures (eg: FusionExpress)
     */

    // Case 4
    if (!_rawMIDs.empty() && !_impl->knowns.empty() && !_impl->splices.empty())
    {
        if (_impl->knowns.size() != _impl->splices.size())
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
            for (auto &j : _impl->splices)
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
         * Constructing differential between normal and fusion genes
         */
        
        for (const auto &i : _impl->normToFus)
        {
            SpliceChimeric s;
            
            const auto normal = match(i.first);
            const auto fusion = match(i.second);

            s.normal = normal->abund(Mix_1);
            s.fusion = fusion->abund(Mix_1);
            
            _impl->spliceChim[i.first]  = s;
            _impl->spliceChim[i.second] = s;
        }
    }

    // Case 5
    else if (!_rawMIDs.empty() && !_impl->splices.empty())
    {
        merge(_rawMIDs, getKeys(_impl->splices));
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
    static BinID createBinID(const GeneID &cID, const GeneID &gID, const IsoformID &iID, const Locus &l)
    {
        return (boost::format("%1%_%2%_%3%_%4%_%5%") % cID
                                                     % gID
                                                     % iID
                                                     % l.start
                                                     % l.end).str();
    }

    struct RawData
    {
        std::map<SequinID, GeneID>                rawMapper;
        std::map<SequinID, std::vector<ExonData>> exonsByTrans;
    };

    struct Data
    {
        // Number of bases for all the reference exons
        Base exonBase;
        
        std::map<GeneID, Locus> _genes;
        std::map<GeneID, GeneData> genes;
        
        std::vector<ExonData>   mergedExons;
        std::vector<ExonData>   sortedExons;
        std::vector<IntronData> sortedIntrons;

        // Only valid for synthetic
        std::map<GeneID, SequinData *> gene2Seqs;
    };
    
    void addRef(const ChrID &cID, const IsoformID &iID, const GeneID &gID, const Locus &l, RawData &raw)
    {
        if (cID != ChrT && !refChrID.empty() && cID != refChrID)
        {
            throw "Multiple chromosomes in the genome not supported";
        }
        
        if (cID != ChrT)
        {
            refChrID = cID;
        }
        
        const auto exon = ExonData(cID, iID, gID, l);

        raw.exonsByTrans[iID].push_back(exon);
        raw.rawMapper[iID] = gID;
    }

    // Eg: chr21...
    ChrID refChrID;
    
    RawData cRaw;

    /*
     * Validated data and resources
     */

    std::map<ChrID, Data> data;
};

TransRef::TransRef() : _impl(new TransRefImpl()) {}

ChrID TransRef::endoID() const
{
    return _impl->refChrID;
}

Base TransRef::exonBase(const ChrID &cID) const
{
    return _impl->data[cID].exonBase;
}

std::vector<GeneID> TransRef::geneIDs(const ChrID &cID) const
{
    std::vector<GeneID> gIDs;

    for (const auto &i : _impl->data.at(cID).genes)
    {
        gIDs.push_back(i.first);
    }

    return gIDs;
}

Limit TransRef::absoluteGene(const SequinHist &hist) const
{
    return Reference<TransData, DefaultStats>::absolute(hist, [&](const GeneID &id)
    {
        return findGene(ChrT, id);
    });
}

void TransRef::addGene(const ChrID &cID, const GeneID &gID, const Locus &l)
{
    // Synthetic is not supported for now...
    assert(cID != ChrT);
    
    if (cID != ChrT)
    {
        _impl->data[cID]._genes[gID] = l;
    }
}

void TransRef::addExon(const ChrID &cID, const GeneID &gID, const IsoformID &iID, const Locus &l)
{
    if (cID == ChrT)
    {
        _impl->addRef(cID, iID, gID, l, _impl->cRaw);
    }
    else
    {
        _impl->data[cID].sortedExons.push_back(ExonData(cID, iID, gID, l));
    }
}

/*
 * ------------------------- Accessors for TransQuin -------------------------
 */

Counts TransRef::countExons(const ChrID &cID) const
{
    return _impl->data.at(cID).sortedExons.size();
}

Counts TransRef::countMerged(const ChrID &cID) const
{
    return _impl->data.at(cID).mergedExons.size();
}

Counts TransRef::countIntrons(const ChrID &cID) const
{
    return _impl->data.at(cID).sortedIntrons.size();
}

const TransRef::GeneData * TransRef::findGene(const ChrID &cID, const GeneID &id) const
{
    return _impl->data.at(cID).genes.count(id) ? &(_impl->data.at(cID).genes.at(id)) : nullptr;
}

template <typename Iter> const typename Iter::mapped_type *findMap(const Iter &x, const Locus &l, MatchRule m)
{
    for (const auto &i : x)
    {
        if ((m == Exact && i.second.l() == l) || (m == Contains && i.second.l().contains(l)))
        {
            return &i.second;
        }
    }
    
    return nullptr;
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

const TransRef::GeneData * TransRef::findGene(const ChrID &cID, const Locus &l, MatchRule m) const
{
    return findMap(_impl->data.at(cID).genes, l, m);
}

const TransRef::ExonData * TransRef::findExon(const ChrID &cID, const Locus &l, MatchRule m) const
{
    return findList(_impl->data.at(cID).sortedExons, l, m);
}

const TransRef::IntronData * TransRef::findIntron(const ChrID &cID, const Locus &l, MatchRule m) const
{
    return findList(_impl->data.at(cID).sortedIntrons, l, m);
}

Intervals<TransRef::ExonInterval> TransRef::exonInters(const ChrID &cID) const
{
    Intervals<ExonInterval> inters;
    
    for (const auto &i : _impl->data.at(cID).sortedExons)
    {
        inters.add(ExonInterval(i.gID, i.iID, TransRefImpl::createBinID(i.cID, i.gID, i.iID, i.l), i.l));
    }
    
    inters.build();
    return inters;
}

Intervals<TransRef::IntronInterval> TransRef::intronInters(const ChrID &cID) const
{
    Intervals<IntronInterval> inters;
    
    for (const auto &i : _impl->data.at(cID).sortedIntrons)
    {
        inters.add(IntronInterval(i.gID, i.iID, TransRefImpl::createBinID(i.cID, i.gID, i.iID, i.l), i.l));
    }

    inters.build();
    return inters;
}

SequinHist TransRef::geneHist(const ChrID &cID) const
{
    if (cID == ChrT)
    {
        return createHist(_impl->data.at(cID).genes);
    }
    else
    {
        return createHist(_impl->data.at(cID)._genes);
    }
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
        data.gID = !_impl->cRaw.rawMapper.empty() ? _impl->cRaw.rawMapper.at(data.id) : "";

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
    
    // Only a single reference chromosome is supported
    assert(_impl->refChrID.size() <= 1);
    
    if (!_impl->data.empty())
    {
        _impl->refChrID = (*_impl->data.begin()).first;
    }
}

template <typename T> void createTrans(const ChrID &cID, T &t)
{
    /*
     * Generate the appropriate structure for analysis
     *
     *   1. Sort the exons
     *   2. Use the sorted exons to generate sorted introns
     *   3. Count the number of non-overlapping bases for the exons
     */

    assert(!t.sortedExons.empty());

    // 1. Sort the exons
    std::sort(t.sortedExons.begin(), t.sortedExons.end(), [](const TransRef::ExonData &x, const TransRef::ExonData &y)
    {
        return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
    });
    
    /*
     * 2. Generate a list of sorted introns, only possible once the exons are sorted.
     */
    
    std::map<SequinID, std::vector<const TransRef::ExonData *>> sorted;
    
    for (const auto &i : t.sortedExons)
    {
        sorted[i.iID].push_back(&i);
    }
    
    for (const auto &i : sorted)
    {
        for (auto j = 1; j < i.second.size(); j++)
        {
            const auto &x = i.second[j-1];
            const auto &y = i.second[j];
            
            TransRef::IntronData d;
            
            d.gID = x->gID;
            d.iID = x->iID;
            d.cID = x->cID;
            d.l   = Locus(x->l.end + 1, y->l.start - 1);
            
            t.sortedIntrons.push_back(d);
        }
    }
    
    // Sort the introns
    std::sort(t.sortedIntrons.begin(), t.sortedIntrons.end(), [](const TransRef::IntronData &x,
                                                                 const TransRef::IntronData &y)
    {
        return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
    });
    
    assert(!t.sortedIntrons.empty());
    
    // 3. Count number of non-overlapping bases for all exons
    t.exonBase = countLocus(t.mergedExons = Locus::merge<TransRef::ExonData, TransRef::ExonData>(t.sortedExons));

    assert(t.exonBase);
}

void TransRef::validate()
{
    const auto iIDs = getKeys(_impl->cRaw.exonsByTrans);
    
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
     * This is only for the synthetic chromosome. Add the validated exons and construct the structures
     * for the genes.
     */
    
    for (const auto &i : _impl->cRaw.exonsByTrans)
    {
        if (_data.count(i.first))
        {
            for (const auto &j : i.second)
            {
                _impl->data[ChrT].sortedExons.push_back(j);
            }
            
            _data[i.first].l = Locus::expand(i.second, [&](const ExonData &f)
            {
                return true;
            });
    
            _impl->data[ChrT].gene2Seqs[_data[i.first].gID] = &_data[i.first];

            _impl->data[ChrT].genes[_data[i.first].gID].id = _data[i.first].gID;         // TODO: ...
            _impl->data[ChrT].genes[_data[i.first].gID].seqs.push_back(&_data[i.first]); // TODO: ...
        }
    }

    /*
     * Create structure for each chromosome
     */
    
    for (const auto &i : _impl->data)
    {
        createTrans(i.first, _impl->data[i.first]);
    }
    
    //assert(_impl->data.count(ChrT));
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
     * Data structure for the genome
     */
    
    ChrID genoID;
    
    // Genomic intervals (eg: chr21)
    Intervals<> inters;

    /*
     * Data structure for the standards
     */
    
    // VarQuin standards
    std::map<SequinID, Locus> stands;

    /*
     * Data structure for the variants
     */
    
    // VarQuin variants (VCF file)
    std::set<Variant> vars;

    /*
     * Data structure for bases (eg: D1_1)
     */

    // Mixture for bases
    std::map<SequinID, Base> baseMix;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

Locus VarRef::matchStand(const SequinID &id) const
{
    return _impl->stands.at(id);
}

Proportion VarRef::matchAlleleFreq(const SequinID &id) const
{
    auto id_ = id;
    
    // Just in case, for example, D_1_1_V
    boost::replace_all(id_, "_V", "_R");

    const auto &p = _impl->data.at(Mix_1).at(id_);
    const auto &r = p.r;
    const auto &v = p.v;

    return v->abund / (r->abund + v->abund);
}

Fold VarRef::matchFold(const SequinID &id) const
{
    const auto &p = _impl->data.at(Mix_1).at(id);
    return round(p.r->abund / p.v->abund);
}

void VarRef::addVar(const Variant &v)
{
    _impl->vars.insert(v);
}

void VarRef::addRInterval(const ChrID &id, const Interval &i)
{
    if (!_impl->genoID.empty() && _impl->genoID != id)
    {
        throw std::runtime_error("Multi chromosomes is not supported. Only a single chromosome can be used.");
    }

    _impl->genoID = id;
    _impl->inters.add(i);
}

void VarRef::addStand(const SequinID &id, const Locus &l)
{
    assert(l.length());
    assert(!_impl->stands.count(id));

    // We're only interested in the position of the sequin
    _impl->stands[id] = l;
}

Counts VarRef::countInters() const
{
    return _impl->inters.size();
}

Counts VarRef::countSeqs() const
{
    return data().size();
}

Counts VarRef::countIndels() const
{
    return std::count_if(_impl->vars.begin(), _impl->vars.end(), [&](const Variant &v)
    {
        return v.type() == Insertion || v.type() == Deletion;
    });
}

Counts VarRef::countSNPs() const
{
    return std::count_if(_impl->vars.begin(), _impl->vars.end(), [&](const Variant &v)
    {
        return v.type() == SNP;
    });
}

Counts VarRef::countVars() const
{
    return _impl->vars.size();
}

void VarRef::validate()
{
    /*
     * Rules:
     *
     *   1: Annotation (eg: VarDiscover)
     *   2: Annotation & mixture (eg: VarAlign)
     *   3: Variants & mixture (eg: VarAllele)
     */
    
    // Rule: 2 and 3
    if (!_rawMIDs.empty())
    {
        merge(_rawMIDs);
    }
    
    // Rule: 1
    else if (!_impl->stands.empty())
    {
        std::set<SequinID> ids;
        
        for (const auto &i : _impl->stands)
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
     * Constructing for the standards
     */
    
    for (const auto &i : _impl->stands)
    {
        _data.at(i.first).l = i.second;
    }

    /*
     * Constructing the genomic intervals (eg: chr21)
     */

    if (_impl->inters.size())
    {
        _impl->inters.build();
    }

    /*
     * Merging the reference and variant sequins. Mixture might not be defined.
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

const VarRef::Base * VarRef::findBase(const SequinID &id, Mixture mix) const
{
    return _impl->baseMix.count(id) ? &(_impl->baseMix.at(id)) : nullptr;
}

Limit VarRef::absoluteBase(const SequinHist &hist, Mixture mix) const
{
    return Reference<SequinData, DefaultStats>::absolute(hist, [&](const SequinID &id)
    {
        return findBase(id);
    }, mix);
}

SequinHist VarRef::baseHist() const
{
    SequinHist hist;
    
    for (const auto &i : _data)
    {
        hist[baseID(i.first)];
    }
    
    return hist;
}

const Intervals<> VarRef::genoInters() const
{
    return _impl->inters;    
}

ChrID VarRef::endoID() const
{
    return _impl->genoID;
}

GenomeHist VarRef::genomeHist() const
{
    GenomeHist hist;
    
    for (const auto &i : _impl->inters.data())
    {
        hist[i.second.id()];
    }

    return hist;
}

Interval * VarRef::findEndo(const Locus &l) const
{
    return _impl->inters.contains(l);
}

const Variant * VarRef::findVar(const SequinID &id) const
{
    // Eg: D_1_1
    auto x = id;
    
    if (boost::algorithm::ends_with(x, "_V"))
    {
        x = x.substr(0, x.length() - 2) + "_R";
    }
    else if (!boost::algorithm::ends_with(x, "_R"))
    {
        x = x + "_R";
    }

    for (const auto &i : _impl->vars)
    {
        if (i.id == x)
        {
            return &i;
        }
    }
    
    return nullptr;
}

const Variant * VarRef::findVar(const Locus &l, MatchRule match) const
{
    if (match != Exact && match != Contains)
    {
        throw std::runtime_error("Only Exact and Contains are supported");
    }

    switch (match)
    {
        case Exact:
        {
            for (const auto &i : _impl->vars)
            {
                if (i.l.start == l.start)
                {
                    return &i;
                }
            }

            break;
        }
            
        case Contains:
        {
            typedef std::pair<SequinID, Locus> StandardPair;
            
            const auto iter = std::find_if(_impl->stands.begin(), _impl->stands.end(), [&](const StandardPair &p)
            {
                return p.second.contains(l);
            });
            
            if (iter != _impl->stands.end())
            {
                return findVar(iter->first);
            }

            break;
        }

        default : { break; }
    }
    
    return nullptr;
}