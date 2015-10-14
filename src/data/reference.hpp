#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include <map>
#include "data/variation.hpp"
#include "stats/sensitivity.hpp"

namespace Anaquin
{
    enum Mixture
    {
        Mix_1,
        Mix_2,
        Mix_3,
        Mix_4,
    };

    struct SequinStats
    {
        // Empty Implementation
    };

    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        // Return the abundance for this sequin specified by the mixture        
        inline Concentration abund(Mixture m, bool norm = false) const
        {
            return mixes.at(m) / (norm ? length : 1);
        }

        SequinID id;

        // Length of the sequin
        Base length;

        // Spiked-in concentration (not available if no mixture provided)
        std::map<Mixture, Concentration> mixes;

        Locus l;
    };

    typedef std::map<SequinID, Counts> SequinHist;

    /*
     * Different rules how two positions can be compared
     */

    enum MatchRule
    {
        Exact,
        Overlap,
        StartOnly,
        Contains,
    };

    template <typename Data = SequinData, typename Stats = SequinStats> class Reference
    {
        public:

            // Add a sequin defined in a mixture file
            void add(const SequinID &id, Base length, Concentration c, Mixture m)
            {
                _mixes[m].insert(MixtureData(id, length, c));
                _rawMIDs.insert(id);
            }

            // Return number of sequins in the mixture
            inline std::size_t countMixes() const { return _mixes.size(); }

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
             * Construct a histogram for each validated sequin
             */

            inline SequinHist hist() const
            {
                SequinHist h;
            
                for (const auto &i : _data)
                {
                    h[i.first] = 0;
                }

                return h;
            }

            // Calculate the total size of all sequins in the reference
            inline Base size() const
            {
                Base n;
                
                for (const auto &i : _data)
                {
                    n += i.second.l.size();
                }
                
                assert(n);
                return n;
            }

            virtual void validate() = 0;

            inline Sensitivity limit(const SequinHist &h) const
            {
                return limit(h, [&](const SequinID &id)
                {
                    return this->match(id);
                });
            }

        protected:

            /*
             * Provides a common mechanism to calculate limit of detection given a histogram or a
             * distribution.
             */

            template <typename F> Sensitivity limit
                                (const std::map<std::string, Counts> &h, F f, Mixture m = Mix_1) const
            {
                Sensitivity s;
            
                // The lowest count must be zero because it can't be negative
                s.counts = std::numeric_limits<unsigned>::max();
            
                for (auto i = h.begin(); i != h.end(); i++)
                {
                    const auto counts = i->second;
                
                    /*
                     * Is this sequin detectable? If it's detectable, what about the concentration?
                     * By definition, detection limit is defined as the smallest abundance while
                     * still being detected.
                     */
                
                    if (counts)
                    {
                        const auto &id = i->first;
                        const auto seq = f(id);

                        // Hard to believe a sequin in the histogram is undefined
                        assert(seq);

                        if (counts < s.counts || (counts == s.counts && seq->abund(m) < s.abund))
                        {
                            s.id     = id;
                            s.counts = counts;
                            s.abund  = seq->abund(m);
                        }
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
                MixtureData(const SequinID &id, Base length, Concentration abund)
                    : id(id), length(length), abund(abund) {}
            
                inline bool operator<(const MixtureData &x)  const { return id < x.id;  }
                inline bool operator==(const MixtureData &x) const { return id == x.id; }
            
                SequinID id;
            
                // Length of the sequin
                Base length;
            
                // Amount of spiked-in concentration
                Concentration abund;
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
                            _data.at(j.id).length = j.length;
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
             * Raw data - structure before validated
             */
        
            // Set of IDs defined in the mixture
            std::set<SequinID> _rawMIDs;

            std::map<Mixture, std::set<MixtureData>> _mixes;
    };

    /*
     * -------------------- Ladder Analysis --------------------
     */
    
    class LadderRef : public Reference<SequinData, SequinStats>
    {
        public:
        
            typedef std::string JoinID;
            typedef std::string UnjoinID;
        
            typedef std::set<JoinID>   JoinIDs;
            typedef std::set<UnjoinID> UnjoinIDs;

            typedef std::map<JoinID, Counts> JoinHist;

            LadderRef();

            void validate() override;

            // Return a list of sequin IDs at the joined level
            JoinIDs joinIDs() const;

            // Construct a histogram at the joined level
            JoinHist joinHist() const;

            // Calculate the limit of detection at the joined level
            Sensitivity limitJoin(const JoinHist &) const;

            // Return abundance for all segments of a particular conjoined
            void abund(const JoinID &, Concentration &, Concentration &, Concentration &, Concentration &,
                             Mixture) const;

        private:
    
            struct LadderRefImpl;
        
            std::shared_ptr<LadderRefImpl> _impl;
    };
    
    /*
     * -------------------- Metagenomics Analysis --------------------
     */
    
    class MetaRef : public Reference<SequinData, SequinStats>
    {
        public:
            MetaRef();
        
            void validate() override;

            void addStand(const SequinID &, Base l);
        
        private:

            struct MetaRefImpl;

            std::shared_ptr<MetaRefImpl> _impl;
    };
    
    /*
     * -------------------- Fusion Analysis --------------------
     */
    
    class FusionRef : public Reference<SequinData, SequinStats>
    {
        public:
        
            struct FusionPoint
            {
                inline bool operator<(const FusionPoint &x)  const { return id < x.id;  }
                inline bool operator==(const FusionPoint &x) const { return id == x.id; }

                operator const SequinID &() const { return id; }
                
                // Where this fusion belongs
                SequinID id;
            
                // The position of the break-point
                Base l1, l2;
            
                // Orientation for each of the segment
                Strand s1, s2;
            };

            FusionRef();

            // Add a known reference fusion
            void addBreak(const FusionPoint &);

            void validate() override;

            // Return the number of reference fusions
            std::size_t countFusions() const;
        
            const FusionPoint * find(Base x, Base y, Strand o1, Strand o2, double fuzzy) const;

        private:
            struct FusionRefImpl;

            std::shared_ptr<FusionRefImpl> _impl;
    };

    /*
     * -------------------- Variant Analysis --------------------
     */
    
    class VarRef : public Reference<SequinData, SequinStats>
    {
        public:

            typedef std::string PairID;
            typedef std::map<PairID, Counts> PairHist;

            /*
             * In VarQuin, we emulate the homozygous and heterozygous genotype
             */
        
            struct GenotypeData
            {
                PairID id;

                // By definition, it should be same as the reference and variant
                Locus l;
                
                const SequinData *r;
                const SequinData *v;

                Concentration abund(Mixture m) const;
            };

            VarRef();

            // Add a reference for a known variant
            void addVar(const Variation &);

            // Add a new standard
            void addStand(const SequinID &, const Locus &);

            // Return number of validated variants
            std::size_t countVars() const;

            // Return number of validated SNPs
            std::size_t countSNPs() const;
        
            // Return number of validated indels
            std::size_t countIndels() const;
        
            // Return number of validated references, eg: D_1_11_R
            std::size_t countRefGenes() const;

            // Return number of validated variants, eg: D_1_11_V
            std::size_t countVarGens() const;
        
            void validate() override;

            // Return the detection limit at the pair level
            Sensitivity limitPair(const PairHist &) const;
       
            // Return a histogram for all the validated pairs
            PairHist pairHist() const;

            // Find a reference gene that contains the given locus
            const GenotypeData *findGeno(const PairID &) const;
        
            // Find a reference gene that contains the given locus
            const GenotypeData *findGeno(const Locus &, double fuzzy = 0, MatchRule = Contains) const;

            // Find a reference variant given a locus
            const Variation *findVar(const Locus &, double fuzzy = 0, MatchRule = StartOnly) const;

            // Return the proportion of variants for a variant pair
            double alleleFreq(Mixture, const PairID &) const;

        private:

            struct VarRefImpl;
            struct VariantPair;

            std::shared_ptr<VarRefImpl> _impl;
    };
    
    /*
     * -------------------- Transcriptome Analysis --------------------
     */
    
    struct TransData : public SequinData
    {
        GeneID gID;
    };

    class TransRef : public Reference<TransData, SequinStats>
    {
        public:

            struct GeneData
            {
                GeneID id;

                inline Locus l() const
                {
                    Base end   = std::numeric_limits<Base>::min();
                    Base start = std::numeric_limits<Base>::max();
                
                    for (const auto &i : seqs)
                    {
                        end   = std::max(end, i->l.end);
                        start = std::min(start, i->l.start);
                    }

                    return Locus(start, end);
                }

                // Calculate the abundance for the gene (summing up all the isoforms)
                inline Concentration abund(Mixture m) const
                {
                    Concentration n = 0;
                
                    for (const auto &i : seqs)
                    {
                        n += (*i).abund(m);
                    }
                
                    return n;
                }
    
                // Each sequin comprises an isoform
                std::vector<SequinData *> seqs;
            };
        
            struct ExonData
            {
                operator const Locus &() const { return l; }

                inline bool operator<(const ExonData &x)  const { return l < x.l;  }
                inline bool operator==(const ExonData &x) const { return l == x.l; }
            
                inline bool operator!=(const Locus &l) const { return !operator==(l); }
                inline bool operator==(const Locus &l) const { return this->l.start == l.start && this->l.end == l.end; }
            
                inline void operator+=(const ExonData &x)
                {
                    l.start = std::min(l.start, x.l.start);
                    l.end   = std::max(l.end, x.l.end);
                }
            
                Locus     l;
                GeneID    gID;
                IsoformID iID;
            };

            struct IntronData
            {
                operator const Locus &() const { return l; }
            
                inline bool operator<(const ExonData &x)  const { return l < x.l;  }
                inline bool operator==(const ExonData &x) const { return l == x.l; }
            
                Locus     l;
                GeneID    gID;
                IsoformID iID;
            };

            typedef std::map<GeneID, Counts> GeneHist;

            TransRef();

            // Return a histogram for all the validated genes
            GeneHist histGene() const;

            // Add a new annoation reference
            void addRef(const IsoformID &iID, const GeneID &gID, const Locus &l);

            void merge(const std::set<SequinID> &mIDs, const std::set<SequinID> &aIDs);

            void validate() override;

            // Calculate the detection limit at the gene level
            Sensitivity limitGene(const GeneHist &) const;

            // Number of non-overlapping bases in all exons
            Base exonBase() const;

            long countMergedExons()   const;
            long countSortedExons()   const;
            long countSortedIntrons() const;
        
            const std::vector<ExonData> & mergedExons() const;

            const GeneData   *findGene  (const GeneID &id)          const;
            const GeneData   *findGene  (const Locus &l, MatchRule) const;
            const ExonData   *findExon  (const Locus &l, MatchRule) const;
            const IntronData *findIntron(const Locus &l, MatchRule) const;

        private:

            struct TransRefImpl;

            std::shared_ptr<TransRefImpl> _impl;        
    };
}

#endif