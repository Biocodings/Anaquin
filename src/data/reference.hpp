#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include "data/hist.hpp"
#include "data/bData.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/minters.hpp"
#include "data/dinters.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    enum class Tool
    {
        Test,
        Help,
        
        RnaAlign,
        RnaAssembly,
        RnaExpress,
        RnaFoldChange,
        RnaSubsample,
        
        VarCopy,
        VarAlign,
        VarDetect,
        VarCancer,
        VarSample,
        VarTrim,
        VarFlip,
        VarKAbund,
        VarStructure,
        VarConjoint,
        
        MetaAlign,
        MetaFoldChange,
        MetaAbund,
        MetaAssembly,
        MetaSubsample
    };
    
    /*
     * Generic template for a sequin. Specalized definitions expected to derive from this class.
     */
    
    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        // Expected concentration
        inline Concent concent(Mixture m = Mix_1, bool norm = false) const
        {
            return mixes.at(m) / (norm ? l.length() : 1);
        }
        
        // Expected fold-change (only valid for two mixtures)
        inline Fold fold() const
        {
            return concent(Mixture::Mix_2) / concent(Mixture::Mix_1);
        }

        SequinID id;

        Locus l;

        // Expected concentration (not available if no mixture provided)
        std::map<Mixture, Concent> mixes;
    };

    /*
     * Different rules how two positions can be compared
     */

    enum MatchRule
    {
        Exact,
        Overlap,
        Contains,
    };
    
    class Ladder;

    struct UserReference
    {
        std::shared_ptr<Ladder> l1, l2;

        // First bed regions (not trimmed)
        std::shared_ptr<BedData> r1;
        
        // Second bed regions (trimmed)
        std::shared_ptr<BedData> r2;
    };
    
    template <typename Data = SequinData> class Reference
    {
        public:

            // All sequins in the reference
            inline std::set<SequinID> seqs() const { return _seqs; }

            // Sequins in ladder 1
            inline std::set<SequinID> l1Seqs() const { return _l1->seqs; }
        
            // Sequins in ladder 2
            inline std::set<SequinID> l2Seqs() const { return _l2->seqs; }
        
            // Concentation in the reference ladder
            inline Concent concent1(const SequinID &id, Mixture m = Mix_1) const
            {
                return _l1->concent(id, m);
            }

            // Concentation in the reference ladder
            inline Concent concent2(const SequinID &id, Mixture m = Mix_1) const
            {
                return _l2->concent(id, m);
            }
        
            // Position in the reference annoation
            inline Locus locus(const SequinID &id) const
            {
                for (const auto &i : *(_r1))
                {
                    if (i.second.r2d.count(id))
                    {
                        return i.second.r2d.at(id).l;
                    }
                }
                
                throw std::runtime_error("Region not found for " + id);
            }
        
            inline Chr2DInters  regs1()  const { return _r1->inters();  }
            inline Chr2DInters  regs2()  const { return _r2->inters();  }
            inline Chr2MInters  mRegs1() const { return _r1->minters(); }
            inline Chr2MInters  mRegs2() const { return _r2->minters(); }
        
            inline Counts nRegs() const { return _r1->count();  }
            inline Counts lRegs() const { return _r1->length(); }
        
            inline MergedIntervals<> mInters(const ChrID &cID) const { return _r1->minters(cID); }
            inline Chr2MInters mInters() const { return _r1->minters(); }

            inline const Data *match(const SequinID &id) const
            {
                return _data.count(id) ? &_data.at(id) : nullptr;
            }

            inline void finalize(Tool x, const UserReference &r)
            {
                validate(x, r);
                
                for (auto &i : _data)
                {
                    if (!i.second.l.length())
                    {
                        throw std::runtime_error("Validation failed. Zero length in data.");
                    }
                }
            }

        protected:

            inline void build(std::shared_ptr<BedData> r1)
            {
                _r1 = r1;
            }

            inline void build(std::shared_ptr<BedData> r1, std::shared_ptr<BedData> r2)
            {
                _r1 = r1;
                _r2 = r2;
            }

            inline void build(std::shared_ptr<Ladder> l1)
            {
                _l1   = l1;
                _seqs = l1->seqs;
            }

            inline void build(std::shared_ptr<Ladder> l1, std::shared_ptr<Ladder> l2)
            {
                _l1 = l1;
                _l2 = l2;
            }

            inline void build(std::shared_ptr<Ladder> l1, std::shared_ptr<BedData> r1)
            {
                _l1   = l1;
                _r1   = r1;
                _seqs = l1->seqs;
            }
        
            inline void build(std::shared_ptr<Ladder>  l1,
                              std::shared_ptr<BedData> r1,
                              std::shared_ptr<BedData> r2)
            {
                _l1   = l1;
                _r1   = r1;
                _r2   = r2;
                _seqs = l1->seqs;
            }

            virtual void validate(Tool, const UserReference &) = 0;

            struct MixtureData
            {
                MixtureData(const SequinID &id, Base length, Concent abund)
                        : id(id), length(length), abund(abund) {}

                inline bool operator<(const SequinID &id)  const { return this->id < id;  }
                inline bool operator==(const SequinID &id) const { return this->id == id; }
                
                inline bool operator<(const MixtureData &x)  const { return id < x.id;  }
                inline bool operator==(const MixtureData &x) const { return id == x.id; }

                SequinID id;

                // Length of the sequin
                Base length;

                // Amount of spiked-in abundance
                Concent abund;
            };

            template <typename T> std::vector<SequinID> merge(const std::set<T> &t1, const std::set<T> &t2)
            {
                std::set<SequinID> x, y;
                
                for (const auto &i : t1) { x.insert(static_cast<SequinID>(i)); }
                for (const auto &i : t2) { y.insert(static_cast<SequinID>(i)); }

                assert(!x.empty() && !y.empty());

                std::vector<SequinID> diffs, inters;
            
                /*
                 * Check for any sequin defined in x but not in y
                 */
            
                std::set_difference(x.begin(),
                                    x.end(),
                                    y.begin(),
                                    y.end(),
                                    std::back_inserter(diffs));

                /*
                 * Check for any sequin defined in both sets
                 */
            
                std::set_intersection(x.begin(),
                                      x.end(),
                                      y.begin(),
                                      y.end(),
                                      std::back_inserter(inters));

                /*
                 * Construct a set of validated sequins. A valid sequin is one in which it's
                 * defined in both mixture and annoation.
                 */
            
                std::for_each(inters.begin(), inters.end(), [&](const SequinID &id)
                {
                    auto d = Data();
                              
                    d.id  = id;
                    
                    // Add a new entry for the validated sequin
                    _data[id] = d;

                    A_ASSERT(!d.id.empty());
                });

                /*
                 * Now, we have a list of validated sequins. Use those sequins to merge information.
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
                            _data.at(j.second->id).mixes[mix] = j.second->abund;
                        }
                    }
                }
            
                A_ASSERT(!_data.empty());
                
                return diffs;
            }
        
            template <typename T> std::vector<SequinID> merge(const std::set<T> &x)
            {
                return merge(x, x);
            }

            // Sequins
            std::set<SequinID> _seqs;
        
            // Sequin regions
            std::shared_ptr<BedData> _r1, _r2;
        
            // Sequin ladder
            std::shared_ptr<Ladder> _l1, _l2;

        
        
        
        
        
        
        
            // Validated ladder
            std::map<Mixture, std::map<SequinID, std::shared_ptr<MixtureData>>> _mixes;

            // Sequins
            std::map<SequinID, Data> _data;
    };

    /*
     * -------------------- Metagenomic Reference --------------------
     */
    
    class MetaRef : public Reference<>
    {
        public:
        
            MetaRef();

            void readBed(const Reader &);

            Counts nMicroSyn() const;
            Counts nMicroGen() const;

            Base nBaseSyn() const;
            Base nBaseGen() const;
        
            Chr2MInters mInters() const;

        protected:
        
            void validate(Tool, const UserReference &) override;
        
        private:

            struct MetaRefImpl;
            std::shared_ptr<MetaRefImpl> _impl;
    };

    /*
     * -------------------- Variant Reference --------------------
     */
    
    class VarRef : public Reference<>
    {
        public:

            VarRef();

            void readVRef(const Reader &);
            void readBRef(const Reader &, Base trim = 0);

            Counts nCNV(int)  const;
            Counts nGeno(Genotype)  const;
            Counts nType(Variation) const;
            Counts nContext(SequinVariant::Context) const;

            // Returns all reference variants
            std::set<Variant> vars() const;

            const SequinVariant &findSeqVar(long) const;
        
            const Variant *findVar(const ChrID &, const Locus &) const;

            // Returns the expected allele frequency
            Proportion findAFreq(const SequinID &) const;

        protected:

            void validate(Tool, const UserReference &) override;

        private:

            struct VarRefImpl;
            struct VariantPair;

            std::shared_ptr<VarRefImpl> _impl;
    };
    
    /*
     * -------------------- Transcriptome Referenceb --------------------
     */
    
    struct GeneData;
    struct TransData;
    
    class RnaRef : public Reference<>
    {
        public:

            RnaRef();

            void readRef(const Reader &);

            std::map<ChrID, Hist> histGene() const;
            std::map<ChrID, Hist> histIsof() const;

            Chr2MInters meInters(Strand str) const;
            Chr2MInters ueInters() const;
            Chr2MInters uiInters() const;

            Base countLenSyn() const;
            Base countLenGen() const;

            Chr2MInters mergedExons() const;
            MergedIntervals<> mergedExons(const ChrID &cID) const;

            // Number of sequin genes from mixture
            Counts nGeneSeqs() const;

            Counts countUExon(const ChrID &) const;
            Counts countUExonSyn() const;
            Counts countUExonGen() const;

            Counts countUIntr(const ChrID &) const;
            Counts countUIntrSyn() const;
            Counts countUIntrGen() const;

            Counts nGeneSyn() const;
            Counts nGeneGen() const;

            Counts countTransSyn() const;
            Counts countTransGen() const;
        
            // Concentration at the gene level
            Concent concent(const GeneID &, Mixture m = Mix_1) const;

            // Expected log-fold at the gene level
            LogFold logFoldGene(const GeneID &) const;

            // Expected log-fold at the sequin (isoform) level
            LogFold logFoldSeq(const IsoformID &) const;

            GeneID s2g(const SequinID &) const;
        
            inline std::set<GeneID> geneIDs() const { return getGenes(ChrIS); }
            inline std::set<IsoformID> isoformIDs() const { return getTrans(ChrIS); }

            std::set<GeneID>  getGenes(const ChrID &) const;
            std::set<TransID> getTrans(const ChrID &) const;

            const GeneData  *findGene (const ChrID &, const GeneID &) const;
            const TransData *findTrans(const ChrID &, const GeneID &) const;

        protected:
        
            void validate(Tool, const UserReference &) override;
            void merge(const std::set<SequinID> &, const std::set<SequinID> &);
        
        private:

            struct RnaRefImpl;

            std::shared_ptr<RnaRefImpl> _impl;        
    };
}

#endif
