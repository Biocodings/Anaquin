#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include <set>
#include "data/data.hpp"
#include "data/hist.hpp"
#include "stats/limit.hpp"
#include "data/variant.hpp"
#include "data/intervals.hpp"

namespace Anaquin
{
    enum Mixture
    {
        Mix_1,
        Mix_2,
    };

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

        // Input concentration
        inline Concent concent(Mixture m = Mix_1, bool norm = false) const
        {
            return mixes.at(m) / (norm ? l.length() : 1);
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

            inline Intervals<> inters() const
            {
                Intervals<> inters;

                for (const auto &i : _data)
                {
                    inters.add(Interval(i.first, i.second.l));
                }
                
                return inters;
            }

            // Add a sequin defined in a mixture file
            inline void add(const SequinID &id, Base length, Concent c, Mixture m)
            {
                _mixes[m].insert(MixtureData(id, length, c));
                _rawMIDs.insert(id);
            }

            inline std::vector<SequinID> seqIDs() const
            {
                std::vector<SequinID> x;
                
                for (const auto &i : _data)
                {
                    x.push_back(i.first);
                }
                
                return x;
            }
        
            // Return all validated sequins
            inline const std::map<SequinID, Data> &data() const { return _data; }

            inline const Data *match(const SequinID &id) const
            {
                return _data.count(id) ? &_data.at(id) : nullptr;
            }

            inline const Data *match(const Locus &l, MatchRule m) const
            {
                for (const auto &i : _data)
                {
                    if ((m == Overlap && i.second.l.overlap(l)) || (m == Contains && i.second.l.contains(l)))
                    {
                        return &i.second;
                    }
                }

                return nullptr;
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

            // Calculate the total length of all sequins in the reference
            inline Base size() const
            {
                Base n = 0;
                
                for (const auto &i : _data)
                {
                    n += i.second.l.length();
                }

                assert(n);
                return n;
            }

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
        
            /*
             * Calculate the absolute detection limits. The limit is defined as the detected sequin with the least abundance.
             */

            inline Limit absolute(const SequinHist &hist) const
            {
                return absolute(hist, [&](const SequinID &id)
                {
                    return this->match(id);
                });
            }

        protected:

            virtual void validate() = 0;

            /*
             * Provides a common mechanism to calculate absolute limit of quantification from a histogram of distribution
             */

            template <typename F> Limit absolute(const std::map<std::string, Counts> &h, F f, Mixture mix = Mix_1) const
            {
                Limit s;
            
                // The lowest count must be zero because it can't be negative
                s.counts = std::numeric_limits<unsigned>::max();
            
                for (auto i = h.begin(); i != h.end(); i++)
                {
                    const auto counts = i->second;
                
                    /*
                     * Is this sequin detectable? If it's detectable, what about the concentration?
                     * By definition, detection limit is defined as the smallest abundance while still detected.
                     */
                
                    if (counts)
                    {
                        const auto &id = i->first;
                       /*
                        const auto seq = f(id);

                        // Hard to believe a sequin in the histogram is undefined
                        assert(seq);

                        if (counts < s.counts || (counts == s.counts && seq->concent(mix) < s.abund))
                        {
                            s.id     = id;
                            s.counts = counts;
                            s.abund  = seq->concent(mix);
                        }
                        */
                    }
                }
            
                if (s.counts == std::numeric_limits<unsigned>::max())
                {
                    s.counts = 0;
                }
            
                return s;
            }
        
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

                    assert(!d.id.empty());
                });

                /*
                 * Now, we have a list of validated sequins. Use those sequins to merge information.
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
        
            template <typename T> void merge(const std::set<T> &x)
            {
                return merge(x, x);
            }

            // Validated sequins
            std::map<SequinID, Data> _data;

            // Statistics about the sequins, only valid after validate()
            Stats _stats;

            /*
             * Raw data - before validation
             */
        
            // Set of IDs defined in the mixture
            std::set<SequinID> _rawMIDs;

            std::map<Mixture, std::set<MixtureData>> _mixes;
    };

    /*
     * -------------------- Structural Analysis --------------------
     */
    
    class StructRef : public Reference<SequinData, DefaultStats>
    {
        public:
        
            void addStruct(const SequinID &) const;
        
        protected:
        
            void validate() override;
        
        private:
        
            struct StructRefImpl;

            std::shared_ptr<StructRefImpl> _impl;
    };
    
    /*
     * -------------------- Ladder Analysis --------------------
     */
    
    class LadderRef : public Reference<SequinData, DefaultStats>
    {
        public:
        
            typedef std::string JoinID;
            typedef std::string UnjoinID;
        
            typedef std::set<JoinID>   JoinIDs;
            typedef std::set<UnjoinID> UnjoinIDs;

            typedef std::map<JoinID, Counts> JoinHist;

            LadderRef();

            // Return a list of sequin IDs at the joined level
            JoinIDs joinIDs() const;

            // Construct a histogram at the joined level
            JoinHist joinHist() const;

            // Calculate the limit of detection at the joined level
            Limit limitJoin(const JoinHist &) const;

            // Return abundance for all segments of a particular conjoined
            void concent(const JoinID &, Concent &, Concent &, Concent &, Concent &, Mixture) const;

        protected:

            void validate() override;

        private:
    
            struct LadderRefImpl;
        
            std::shared_ptr<LadderRefImpl> _impl;
    };
    
    /*
     * -------------------- Metagenomics Analysis --------------------
     */
    
    class MetaRef : public Reference<SequinData, DefaultStats>
    {
        public:
            MetaRef();

            void addStand(const SequinID &, const Locus &);

            // Whether the locus is contained in one of the genomes
            const SequinData * contains(const GenomeID &, const Locus &) const;

        protected:

            void validate() override;
        
        private:
        
            struct MetaRefImpl;

            std::shared_ptr<MetaRefImpl> _impl;
    };
    
    /*
     * -------------------- Fusion Analysis --------------------
     */
    
    class FusionRef : public Reference<SequinData, DefaultStats>
    {
        public:

            struct KnownFusion
            {
                inline bool operator<(const KnownFusion &x)  const { return id < x.id;  }
                inline bool operator==(const KnownFusion &x) const { return id == x.id; }

                operator const SequinID &() const { return id; }
                
                // Where this fusion belongs
                SequinID id;
            
                // The position of the breakpoint
                Base l1, l2;
            
                // Orientation for each of the segment
                Strand s1, s2;
            };

            /*
             * Represents the concentration between normal and fusion gene
             */

            struct NormalFusion
            {
                // Concentration for the normal splicing
                Concent normal;
                
                // Concentration for the fusion chimeria
                Concent fusion;

                inline Fold fold() const { return normal / fusion; }
            };

            FusionRef();

            /*
             * Modifier operations
             */
        
            void addFusion(const KnownFusion &);
        
            // Add an intron junction for the normal genes
            void addJunct(const SequinID &, const Locus &);

            // Add fusion or normal genes comprised of the standards
            void addStand(const SequinID &, const Locus &);

            SequinHist normalHist() const;
            SequinHist fusionHist() const;

            Counts countFusion() const;
            Counts countJuncts() const;

            const KnownFusion *findFusion(const SequinID &) const;
            const KnownFusion *findFusion(Base x, Base y, Strand o1, Strand o2, double fuzzy) const;

            // Find a reference junction for the normal genes
            const SequinData *findJunct(const Locus &) const;

            const NormalFusion *findNormFus(const SequinID &) const;

        protected:
        
            void validate() override;

        private:
        
            struct FusionRefImpl;

            std::shared_ptr<FusionRefImpl> _impl;
    };

    /*
     * -------------------- Variant Analysis --------------------
     */
    
    struct Variant;
    
    class VarRef : public Reference<SequinData, DefaultStats>
    {
        public:

            struct Base : public Matched
            {
                // Eg: D_1_1 and D_2_2
                SequinID id;
                
                inline SequinID name() const override
                {
                    return id;
                }

                inline Concent concent(Mixture mix = Mix_1) const
                {
                    return total.at(mix);
                }
                
                std::map<Mixture, Concent> total;
            };

            VarRef();

            /*
             * Modifier functions
             */
        
            // Add a known variant (synthetic and genomic/user)
            void addVar(const Variant &);
        
            // Add a sequin in the standards
            void addStand(const SequinID &, const Locus &);

            // Add a reference interval (eg: chr21)
            void addGInterval(const ChrID &, const Interval &);

            /*
             * Query functions
             */
        
            ChrID genoID() const;

            bool isGenoID(const ChrID &cID) const { return cID == genoID(); }

            const Intervals<> genoInters() const;
        
            // Absolute detection limit at the base level
            Limit absoluteBase(const SequinHist &, Mixture mix = Mix_1) const;

            // Returns number of known variants
            Counts countVars() const;

            // Count SNPs for a chromosome
            Counts countSNP(const ChrID &) const;
        
            // Counts indels for the synthetic chromosomes
            Counts countSNPSync() const;
        
            // Counts indels for the genome
            Counts countSNPGeno() const;
        
            // Counts indels for a chromosome
            Counts countIndel(const ChrID &) const;
        
            // Counts indels for the synthetic chromosomes
            Counts countIndSync() const;
        
            // Counts indels for the genome
            Counts countIndGeno() const;

            inline Counts countSync() const { return countSNPSync() + countIndSync(); }

            // Returns number of sequins
            Counts countSeqs() const;

            // Returns number of genomic intervals
            Counts countInters() const;

            // Eg: D_1_11, D_1_12
            SequinHist baseHist() const;

            // Eg: chr21 intervals
            GenomeHist genomeHist() const;
        
            HashHist varHist() const;
        
            Counts countIntervals(const ChrID &) const;
        
            /*
             * Finding functions
             */
        
            const Variant *hashVar(long key) const;
        
            const Variant *findVar(const SequinID &) const;
            const Variant *findVar(const Locus &, MatchRule = Exact) const;

            Interval *findGeno(const ChrID &, const Locus &, MatchRule) const;
        
            Interval *findGeno(const Locus &) const;
            Interval *findGeno(const ChrID &cID, const Locus &l) const
            {
                return isGenoID(cID) ? findGeno(l) : nullptr;
            }
        
            const Base *findGene(const SequinID &, Mixture mix = Mix_1) const;

            // Returns the expected allele fold-change
            Fold findAFold(const SequinID &) const;

            // Returns the expected allele frequency
            Proportion findAFreq(const SequinID &) const;

        protected:

            void validate() override;

        private:

            struct VarRefImpl;
            struct VariantPair;

            std::shared_ptr<VarRefImpl> _impl;
    };
    
    class Reader;
    
    /*
     * -------------------- Transcriptome Analysis --------------------
     */
    
    struct TransData : public SequinData
    {
        GeneID gID;
    };

    class TransRef : public Reference<TransData, DefaultStats>
    {
        public:

            TransRef();
        
            void readRef(const Reader &);

            // Returns histogram for genes for all chromosomes (synthetic + genome)
            std::map<ChrID, Hist> histGene() const;

            SequinHist geneHist(const ChrID &) const;

            // Intervals for reference exons
            Intervals<> exonInters(const ChrID &) const;
        
            // Intervals for reference introns
            Intervals<> intronInters(const ChrID &) const;
        
            // Absolute detection limit at the gene level
            Limit absoluteGene(const SequinHist &) const;

            Base countLenSyn() const;
            Base countLenGen() const;

            Counts countExonSyn() const;
            Counts countExonGen() const;
        
            Counts countIntrSyn() const;
            Counts countIntrGen() const;

            Counts countGeneSyn() const;
            Counts countGeneGen() const;

            Counts countTransSyn() const;
            Counts countTransGen() const;
        
            // Concentration at the gene level
            Concent concent(const GeneID &, Mixture m = Mix_1) const;
        
            const GeneData *findGene(const ChrID &, const GeneID &)           const;
            const Interval *findGene(const ChrID &, const Locus &, MatchRule) const;

        protected:
        
            void validate() override;

            void merge(const std::set<SequinID> &, const std::set<SequinID> &);
        
        private:

            struct TransRefImpl;

            std::shared_ptr<TransRefImpl> _impl;        
    };
}

#endif