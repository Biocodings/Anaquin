#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include "data/hist.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/minters.hpp"
#include "data/intervals.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct DefaultStats
    {
        // Empty Implementation
    };

    /*
     * Generic template for a sequin. Specalized definitions expected to derive from this class.
     */
    
    struct SequinData : public Matched
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        // Expected concentration
        inline Concent concent(Mixture m = Mix_1, bool norm = false) const
        {
            return mixes.at(m) / (norm ? l.length() : 1);
        }
        
        // Expected differential fold-change (only defined for two mixtures)
        inline Fold fold() const
        {
            return concent(Mixture::Mix_2) / concent(Mixture::Mix_1);
        }

        inline SequinID name() const override
        {
            return id;
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
    
    template <typename Data = SequinData, typename Stats = DefaultStats> class Reference
    {
        public:

            // Add a sequin defined in a mixture file
            inline void add(const SequinID &id, Base length, Concent c, Mixture m)
            {
                _mixes[m][id] = std::shared_ptr<MixtureData>(new MixtureData(id, length, c));
                _rawMIDs.insert(id);
            }

            // Number of sequins defined in mixture
            inline Counts countMix() const { return _mixes.size(); }

            // Number of sequins defined in annotation
            inline Counts countSeqs() const { return _data.size(); }

            // Return all validated sequins
            inline const std::map<SequinID, Data> &data() const { return _data; }

            inline const Data *match(const SequinID &id) const
            {
                return _data.count(id) ? &_data.at(id) : nullptr;
            }

            /*
             * Histogram for distribution
             */

            template <typename T> Hist hist(const T &data) const
            {
                Hist hist;
            
                for (const auto &i : data)
                {
                    hist[i.first] = 0;
                }

                return hist;
            }

            inline SequinHist hist() const { return hist(_data); }

            inline void finalize()
            {
                validate();
                
                for (auto &i : _data)
                {
                    if (!i.second.l.length())
                    {
                        throw std::runtime_error("Validation failed. Zero length in data.");
                    }
                }
            }

        protected:

            virtual void validate() = 0;

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

            /*
             * Provide a common framework for validation. Typically, the sequins can be validated
             * by two set of IDs, for example, mixtute and annotation. This function can also be
             * validated a single set of sequins, simply call merge(x, x).
             */
        
            template <typename T> void merge(const std::set<T> &t1, const std::set<T> &t2)
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
            }
        
            template <typename T> void merge(const std::set<T> &x)
            {
                return merge(x, x);
            }

            inline const MixtureData * findMix(Mixture mix, const SequinID &id) const
            {
                return _mixes.at(mix).count(id) ? _mixes.at(mix).at(id).get() : nullptr;
            }

            // Validated sequins
            std::map<SequinID, Data> _data;

            // Set of IDs defined in the mixture
            std::set<SequinID> _rawMIDs;

            // Data for mixture (if defined)
            std::map<Mixture, std::map<SequinID, std::shared_ptr<MixtureData>>> _mixes;
    };

    /*
     * -------------------- Metagenomic Analysis --------------------
     */
    
    class MetaRef : public Reference<SequinData, DefaultStats>
    {
        public:
        
            MetaRef();

            void readBed(const Reader &);

            Counts nMicroSyn() const;
            Counts nMicroGen() const;

            Base nBaseSyn() const;
            Base nBaseGen() const;
        
            MC2Intervals mInters() const;

        protected:
        
            void validate() override;
        
        private:

            struct MetaRefImpl;
            std::shared_ptr<MetaRefImpl> _impl;
    };

    /*
     * -------------------- Variant Analysis --------------------
     */
    
    class VarRef : public Reference<SequinData, DefaultStats>
    {
        public:

            VarRef();

            void readVRef(const Reader &);
            void readGBRef(const Reader &, Base trim = 0);

            // Number of reference regions
            Counts nRegs() const;

            // Total number of bases for all reference regions
            Base lRegs() const;

            Counts nCNV(int)  const;
            Counts nGeno(Genotype)  const;
            Counts nType(Variation) const;
            Counts nContext(SeqVariant::Context) const;

            Counts countSNP() const;
            Counts countInd() const;
        
            C2Intervals regions(bool trimmed = false) const;

            MC2Intervals mInters() const;
        
            MergedIntervals<> mInters(const ChrID &) const;

            // Returns all reference variants
            std::set<Variant> vars() const;

            const SeqVariant &findSeqVar(long) const;
        
            const Variant *findVar(const ChrID &, long) const;
            const Variant *findVar(const ChrID &, const Locus &) const;

            // Returns the expected allele frequency
            Proportion findAFreq(const SequinID &) const;

        protected:

            void validate() override;

        private:

            struct VarRefImpl;
            struct VariantPair;

            std::shared_ptr<VarRefImpl> _impl;
    };
    
    /*
     * -------------------- Transcriptome Analysis --------------------
     */
    
    struct GeneData;
    struct TransData;
    
    class RnaRef : public Reference<SequinData, DefaultStats>
    {
        public:

            RnaRef();

            void readRef(const Reader &);

            std::map<ChrID, Hist> histGene() const;
            std::map<ChrID, Hist> histIsof() const;

            MC2Intervals meInters(Strand str) const;
            MC2Intervals ueInters() const;
            MC2Intervals uiInters() const;

            Base countLenSyn() const;
            Base countLenGen() const;

            MC2Intervals mergedExons() const;
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
        
            void validate() override;

            void merge(const std::set<SequinID> &, const std::set<SequinID> &);
        
        private:

            struct RnaRefImpl;

            std::shared_ptr<RnaRefImpl> _impl;        
    };
}

#endif
