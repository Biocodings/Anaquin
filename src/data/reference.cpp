#include "data/tokens.hpp"
#include "data/reference.hpp"
#include <iostream>
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
        
        inline Concentration abund(Mixture m) const
        {
            return A->abund(m) + B->abund(m) + C->abund(m) + D->abund(m);
        }
    };

    /*
     * Validated data
     */
    
    std::map<JoinID, JoinData> joined;

    /*
     * Raw data
     */
    
    std::map<JoinID, JoinData> rawJoined;
};

LadderRef::LadderRef() : _impl(new LadderRefImpl()) {}

Limit LadderRef::limitJoin(const JoinHist &h) const
{
    return Reference<SequinData, SequinStats>::limit(h, [&](const JoinID &id)
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
     
        if (typeID == "A") { _impl->joined[joinID].A = &i.second; }
        if (typeID == "B") { _impl->joined[joinID].B = &i.second; }
        if (typeID == "C") { _impl->joined[joinID].C = &i.second; }
        if (typeID == "D") { _impl->joined[joinID].D = &i.second; }
    }

    assert(!_impl->joined.empty());
}

void LadderRef::abund(const LadderRef::JoinID &id, Concentration &a, Concentration &b, Concentration &c,
                                                   Concentration &d, Mixture m) const
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
    /*
     * Validated variables
     */

    // Normal splicing
    std::map<SequinID, Locus> splice;

    // Fusion breaks
    std::set<FusionPoint> breaks;

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
    
    /*
     * Raw variables
     */

    // Standards (fusion genes will have two loci)
    std::map<SequinID, std::vector<Locus>> rawStands;
    
    // Fusion breaks
    std::set<FusionPoint> rawBreaks;
    
    // Normal splicing
    std::map<SequinID, Locus> rawSplices;
};

FusionRef::FusionRef() : _impl(new FusionRefImpl()) {}

Counts FusionRef::countFusion() const { return _impl->breaks.size(); }
Counts FusionRef::countSplice() const { return _impl->splice.size(); }

const FusionRef::SpliceChimeric * FusionRef::findSpliceChim(const SequinID &id) const
{
    return _impl->spliceChim.count(id) ? &_impl->spliceChim.at(id) : nullptr;
}

SequinID FusionRef::normalToFusion(const SequinID &id) const
{
    return _impl->normToFus.at(id);
}

void FusionRef::addBreak(const FusionPoint &f)
{
    _impl->rawBreaks.insert(f);
}

void FusionRef::addSplice(const SequinID &id, const Locus &l)
{
    _impl->rawSplices[id] = l;
}

void FusionRef::addStand(const SequinID &id, const Locus &l)
{
    _impl->rawStands[id].push_back(l);
}

const SequinData *FusionRef::findFusion(const Locus &l) const
{
    assert(!_impl->fusions.empty());
    
    for (auto &i : _impl->fusions)
    {
        if (i.second[0].contains(l) || i.second[1].contains(l))
        {
            return &(_data.at(i.first));
        }
    }
    
    return nullptr;
}

const SequinData *FusionRef::findSplice(const Locus &l) const
{
    assert(!_impl->splice.empty());
    
    for (auto &i : _impl->splice)
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
     *   5: Splicing & Mixtures (eg: FusionNormal)
     */

    // Case 4
    if (!_rawMIDs.empty() && !_impl->rawBreaks.empty() && !_impl->rawSplices.empty())
    {
        if (_impl->rawBreaks.size() != _impl->rawSplices.size())
        {
            throw std::runtime_error("Number of fusions not equal to splicing. Please check and try again.");
        }

        merge(_rawMIDs);
        
        _impl->splice = _impl->rawSplices;
        
        /*
         * Constructing a mapping between splicing and fusion. While we can't assume the orders, we can
         * assume like: "FG1_12_P2" and "NG1_12_P2".
         */
        
        auto f = [&](const SequinID &x, const SequinID &y)
        {
            return x.substr(2, x.size()-1) == y.substr(2, y.size()-1);
        };
        
        for (auto &i: _impl->rawBreaks)
        {
            for (auto &j : _impl->rawSplices)
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
        
        if (_impl->fusToNorm.size() != _impl->rawBreaks.size())
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
    else if (!_rawMIDs.empty() && !_impl->rawSplices.empty())
    {
        merge(_rawMIDs, getKeys(_impl->rawSplices));
        _impl->splice = _impl->rawSplices;
    }

    // Case 3
    else if (!_rawMIDs.empty() && !_impl->rawStands.empty())
    {
        merge(_rawMIDs);

        for (const auto &i : _impl->rawStands)
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
    else if (!_impl->rawBreaks.empty())
    {
        merge(_impl->rawBreaks);
    }

    else
    {
        throw std::runtime_error("Unknown validation");
    }

    /*
     * Copy fusion points (no harm if none provided)
     */
    
    for (const auto &i : _impl->rawBreaks)
    {
        if (_data.count(i.id))
        {
            _impl->breaks.insert(i);
        }
    }
}

inline bool compare(Base x, Base y, Base fuzzy = 0.0)
{
    return std::abs(x - y) <= fuzzy;
}

const Anaquin::FusionRef::FusionPoint * FusionRef::find(Base x, Base y, Strand o1, Strand o2, double fuzzy) const
{
    for (const auto &f : _impl->breaks)
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
    /*
     * Validated variables
     */
    
    // Number of bases for all the reference exons
    Base exonBase;

    std::map<GeneID, GeneData> genes;
    std::vector<ExonData>      mergedExons;
    std::vector<ExonData>      sortedExons;
    std::vector<IntronData>    sortedIntrons;

    /*
     * Raw variables
     */
    
    std::set<GeneID>                           rawGIDs;
    std::set<IsoformID>                        rawIIDs;
    std::map<SequinID, GeneID>                 rawMapper;
    std::map<GeneID, std::vector<ExonData>>    exonsByGenes;
    std::map<IsoformID, std::vector<ExonData>> exonsByTrans;
};

TransRef::TransRef() : _impl(new TransRefImpl()) {}

Base TransRef::exonBase() const { return _impl->exonBase; }

Limit TransRef::limitGene(const GeneHist &h) const
{
    return Reference<TransData, SequinStats>::limit(h, [&](const GeneID &id)
    {
        return findGene(id);
    });
}

void TransRef::addRef(const IsoformID &iID, const GeneID &gID, const Locus &l)
{
    TransRef::ExonData e;
    
    e.l   = l;
    e.iID = iID;
    e.gID = gID;
    
    _impl->exonsByGenes[gID].push_back(e);
    _impl->exonsByTrans[iID].push_back(e);
    
    _impl->rawIIDs.insert(iID);
    _impl->rawGIDs.insert(gID);
    _impl->rawMapper[iID] = gID;
}

long TransRef::countSortedExons()   const { return _impl->sortedExons.size();   }
long TransRef::countMergedExons()   const { return _impl->mergedExons.size();   }
long TransRef::countSortedIntrons() const { return _impl->sortedIntrons.size(); }

const std::vector<TransRef::ExonData> & TransRef::mergedExons()   const { return _impl->mergedExons; }

const TransRef::GeneData * TransRef::findGene(const GeneID &id) const
{
    return _impl->genes.count(id) ? &(_impl->genes.at(id)) : nullptr;
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

const TransRef::GeneData * TransRef::findGene(const Locus &l, MatchRule m) const
{
    return findMap(_impl->genes, l, m);
}

const TransRef::ExonData * TransRef::findExon(const Locus &l, MatchRule m) const
{
    return findList(_impl->sortedExons, l, m);
}

const TransRef::IntronData * TransRef::findIntron(const Locus &l, MatchRule m) const
{
    return findList(_impl->sortedIntrons, l, m);
}

TransRef::GeneHist TransRef::histGene() const
{
    GeneHist h;
    
    for (const auto &i : _impl->genes)
    {
        h[i.first] = 0;
    }
    
    return h;
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
        auto d = TransData();
        
        d.id  = id;
        d.gID = _impl->rawMapper.at(d.id);
        
        // Add a new entry for the validated sequin
        _data[id] = d;
        
        assert(!d.id.empty() && !d.gID.empty());
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
                _data.at(j.id).length = j.length;
                _data.at(j.id).mixes[mix] = j.abund;
            }
        }
    }
    
    assert(!_data.empty());
    
    /*
     * Compute the locus for each sequin
     */
    
    for (const auto &i : _impl->exonsByTrans)
    {
        if (_data.count(i.first))
        {
            _data[i.first].l = Locus::expand(i.second, [&](const ExonData &f)
            {
                return true;
            });

            _impl->genes[_data[i.first].gID].id = _data[i.first].gID;
            _impl->genes[_data[i.first].gID].seqs.push_back(&_data[i.first]);
        }
    }
}

void TransRef::validate()
{
    if (_impl->exonsByGenes.empty())
    {
        throw std::runtime_error("There is no synthetic chromosome in the annotation file. Anaquin is unable to proceed unless a valid annotation is given. Please check your file and try again.");
    }
    
    merge(_rawMIDs, _impl->rawIIDs);
    
    /*
     * Filter out only those validated exons
     */
    
    for (const auto &i : _impl->exonsByTrans)
    {
        if (_data.count(i.first))
        {
            for (const auto &j : i.second)
            {
                _impl->sortedExons.push_back(j);
            }
        }
    }
    
    /*
     * Generate a list of sorted exons
     */
    
    assert(!_impl->sortedExons.empty());
    std::sort(_impl->sortedExons.begin(), _impl->sortedExons.end(), [](const ExonData &x, const ExonData &y)
    {
        return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
    });

    /*
     * Generate a list of sorted introns, only possible once the exons are sorted.
     */
    
    for (auto i = 0; i < _impl->sortedExons.size(); i++)
    {
        if (i && _impl->sortedExons[i].gID == _impl->sortedExons[i-1].gID)
        {
            IntronData d;
            
            d.gID = _impl->sortedExons[i].gID;
            d.iID = _impl->sortedExons[i].iID;
            d.l   = Locus(_impl->sortedExons[i-1].l.end + 1, _impl->sortedExons[i].l.start - 1);
            
            _impl->sortedIntrons.push_back(d);
        }
    }
    
    assert(!_impl->sortedIntrons.empty());
    std::sort(_impl->sortedIntrons.begin(), _impl->sortedIntrons.end(),
                [](const IntronData &x, const IntronData &y)
                {
                    return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
                });

    // Count number of non-overlapping bases for all exons
    _impl->exonBase = countLocus(_impl->mergedExons = Locus::merge<ExonData, ExonData>(_impl->sortedExons));
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
    /*
     * Validated variables
     */

    // Validated variants
    std::set<Variation> vars;

    std::set<SequinID> varIDs;
    
    // Validated genotypes (references + variants)
    std::map<GenoID, GenotypeData> genos;

    std::map<Mixture, std::map<GenoID, VariantPair>> pairs;

    // Reference intervals (eg: chr21)
    std::map<ChromoID, Intervals> inters;

    /*
     * Raw variables
     */
    
    std::set<Variation> rawVars;
    std::set<SequinID>  rawVarIDs;

    // Reference intervals (eg: chr21)
    std::map<ChromoID, Intervals> rawInters;
    
    // Locus of the sequin
    std::map<SequinID, Locus> rawSeqsByID;
};

VarRef::VarRef() : _impl(new VarRefImpl()) {}

double VarRef::alleleFreq(Mixture m, const GenoID &bID) const
{
    const auto &p = _impl->pairs.at(m).at(bID);
    const auto &r = p.r;
    const auto &v = p.v;

    // Abundance ratio of reference to variant DNA standard
    return v->abund / (r->abund + v->abund);
}

void VarRef::addVar(const Variation &v)
{
    assert(!v.bID.empty());
    _impl->rawVars.insert(v);
    _impl->rawVarIDs.insert(v.id);
}

void VarRef::addInterval(const ChromoID &id, const Interval &i)
{
    _impl->rawInters[id].add(i);
}

void VarRef::addStand(const SequinID &id, const Locus &l)
{
    assert(l.length());
    assert(!_impl->rawSeqsByID.count(id));

    // We're only interested in the position of the sequin
    _impl->rawSeqsByID[id] = l;
}

std::size_t VarRef::countRefGenes() const
{
    return static_cast<std::size_t>(0.5 * data().size());
}

std::size_t VarRef::countVarGens() const
{
    return countRefGenes();
}

std::size_t VarRef::countIndels() const
{
    return std::count_if(_impl->vars.begin(), _impl->vars.end(), [&](const Variation &v)
    {
        return v.type == Insertion || v.type == Deletion;
    });
}

std::size_t VarRef::countSNPs() const
{
    return std::count_if(_impl->vars.begin(), _impl->vars.end(), [&](const Variation &v)
    {
        return v.type == SNP;
    });
}

std::size_t VarRef::countVars() const
{
    return _impl->vars.size();
}

VarRef::GenoHist VarRef::genoHist() const
{
    GenoHist h;
    
    for (const auto &i : _impl->genos)
    {
        h[i.first] = 0;
    }
    
    return h;
}

const VarRef::GenotypeData * VarRef::findGeno(const GenoID &id) const
{
    return _impl->genos.count(id) ? &(_impl->genos.at(id)) : nullptr;
}

const VarRef::GenotypeData * VarRef::findGeno(const Locus &l, double fuzzy, MatchRule m) const
{
    for (const auto &i : _impl->genos)
    {
        if ((m == Overlap && i.second.l.overlap(l)) || (m == Contains && i.second.l.contains(l)))
        {
            return &i.second;
        }
    }
    
    return nullptr;
}

Limit VarRef::limitGeno(const GenoHist &h) const
{
    return Reference<SequinData, SequinStats>::limit(h, [&](const VarRef::GenoID &id)
    {
        return findGeno(id);
    });
}

Concentration VarRef::GenotypeData::abund(Mixture m) const
{
    return r->abund(m) + v->abund(m);
}

void VarRef::validate()
{
    _impl->vars   = _impl->rawVars;
    _impl->varIDs = _impl->rawVarIDs;
    _impl->inters = _impl->rawInters;

    /*
     * Validation rules:
     *
     *   1: Annotation (eg: VarCoverage)
     *   2: Annotation & mixture (eg: VarAlign)
     *   3: Variants (eg: VarDiscover)
     *   4: Variants & mixture (eg: VarAllele)
     */
    
    // Rule: 2 and 4
    if (!_rawMIDs.empty())
    {
        // Validate by mixture
        merge(_rawMIDs);
    }
    
    // Rule: 3
    else if (!_impl->varIDs.empty())
    {
        // Validate by variants
        merge(_impl->varIDs);
    }
    
    // Rule: 1
    else if (!_impl->rawSeqsByID.empty())
    {
        std::set<SequinID> ids;
        
        for (const auto &i : _impl->rawSeqsByID)
        {
            ids.insert(i.first);
        }

        // Validate by annotation
        merge(ids);
    }
    else
    {
        throw std::runtime_error("Failed to validate for VarQuin");
    }
    
    /*
     * Construct data-structure for the standards. The purpose is to combine reference and variant
     * sequins. Although it can also be done for the variants, but it's probably not a good idea.
     *
     * Consider:
     *
     *     chrT    641707  641708  D_1_10_R_G/A    0       +
     *     chrT    641714  641715  D_1_10_R_G/A    0       +
     *
     * There is no information directly for D_1_1_V. Therefore, we'll only do it for the standards.
     */
    
    for (const auto &i : _impl->rawSeqsByID)
    {
        if (_data.count(i.first))
        {
            _data[i.first].l = i.second;

            const auto pairID = i.first.substr(0, i.first.size() - 2);

            // TODO: Do we need this?
            _data[i.first].length = i.second.length();
            
            _impl->genos[pairID].l  = i.second;
            _impl->genos[pairID].id = pairID;
            _impl->genos[pairID].r  = &(_data.at(pairID + "_R"));
            _impl->genos[pairID].v  = &(_data.at(pairID + "_V"));
        }
    }
    
    /*
     * Construct data structure for homozygous/heterozygous
     */

    std::vector<std::string> toks;
    
    for (const auto &i : _mixes)
    {
        for (const auto &j : _mixes.at(i.first))
        {
            Tokens::split(j.id, "_", toks);

            // It has be the reference or variant...
            assert(toks[3] == "R" || toks[3] == "V");
            
            // Eg: D_1_10
            const auto baseID = toks[0] + "_" + toks[1] + "_" + toks[2];

            if (toks[3] == "R")
            {
                _impl->pairs[i.first][baseID].r = &j;
            }
            else
            {
                _impl->pairs[i.first][baseID].v = &j;
            }
        }
    }
}

const Interval * VarRef::findQuery(const ChromoID &chr, const Locus &l) const
{
    if (!_impl->inters.count(chr))
    {
        throw std::runtime_error("Failed to find " + chr + " in the reference intervals");
    }

    return _impl->inters.at(chr).contains(l);
}

const Variation * VarRef::findVar(const Locus &l, double fuzzy, MatchRule match) const
{
    for (const auto &i : _impl->vars)
    {
        if (match == StartOnly && i.l.start == l.start)
        {
            return &i;
        }
    }

    return nullptr;
}