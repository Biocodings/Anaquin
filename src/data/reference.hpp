#ifndef GI_REFERENCE_HPP
#define GI_REFERENCE_HPP

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
    };

    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        // Return the normalized abundance
        inline Concentration abund(Mixture m) const
        {
            return mixes.at(m) / length;
        }

        SequinID id;

        // Length of the sequin
        Base length;

        // Amount of spiked-in concentration
        std::map<Mixture, Concentration> mixes;

        Locus l;
    };
    
    typedef std::map<SequinID, Counts> SequinHist;

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

            inline const Data *seq(const SequinID &id) const
            {
                return _data.count(id) ? &_data.at(id) : nullptr;
            }

            inline const Data *seq(const Locus &l) const
            {
                for (const auto &i : _data)
                {
                    if (i.second.l.overlap(l))
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

            virtual void validate() = 0;

            inline Sensitivity limit(const SequinHist &h) const
            {
                return limit(h, [&](const SequinID &id)
                {
                    return this->seq(id);
                });
            }

        protected:

            /*
             * Provides a common mechanism to calculate limit of detection given a histogram or a distribution
             */

            template <typename F> Sensitivity limit(const std::map<std::string, Counts> &h, F f) const
            {
                Sensitivity s;
            
                // The lowest count must be zero because it can't be negative
                s.counts = std::numeric_limits<unsigned>::max();
            
                for (auto iter = h.begin(); iter != h.end(); iter++)
                {
                    const auto counts = iter->second;
                
                    /*
                     * Is this sequin detectable? If it's detectable, what about the concentration?
                     * By definition, detection limit is defined as the smallest abundance while
                     * still being detected.
                     */
                
                    if (counts)
                    {
                        const auto &id = iter->first;
                        const auto seq = f(id);
                    
                        // Hard to believe a sequin in the histogram is undefined
                        assert(seq);
                    
                        if (counts < s.counts || (counts == s.counts && seq->abund(Mix_1) < s.abund))
                        {
                            s.id     = id;
                            s.counts = counts;
                            s.abund  = seq->abund(Mix_1);
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
             * Most sequins can be validated by merging two sets of IDs, typically, a set of
             * IDs defined in the mixture and a set of IDs defined for the annotation. This
             * function provides a common framework for merging the two sets.
             */

            inline void merge(const std::set<SequinID> &x, const std::set<SequinID> &y)
            {
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
    
    struct LadderData : public SequinData
    {
    };

    class LadderRef : public Reference<LadderData, SequinStats>
    {
        public:
        
            inline void validate() override
            {
            
            }
    };
    
    /*
     * -------------------- Metagenomics Analysis --------------------
     */
    
    struct MetaData : public SequinData
    {
    };

    class MetaRef : public Reference<MetaData, SequinStats>
    {
        public:
        
            inline void validate() override
            {
            
            }
    };
    
    /*
     * -------------------- Fusion Analysis --------------------
     */
    
    struct FusionData : public SequinData
    {
        
    };

    class FusionRef : public Reference<FusionData, SequinStats>
    {
        public:
        
            inline void validate() override
            {
            
            }
    };

    /*
     * -------------------- Variant Analysis --------------------
     */
    
    class VarRef : public Reference<SequinData, SequinStats>
    {
        public:

            enum Matching
            {
                StartOnly,
            };

            VarRef();

            // Add a reference for a known variant
            void addVar(const Variation &v);

            // Return the number of validated known variants
            std::size_t countVars() const;

            void validate() override;

            const Variation *findVar(const Locus &l, double fuzzy = 0, Matching = StartOnly) const;

            double alleleFreq(Mixture m, const BaseID &) const;
        
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

            enum Matching
            {
                Exact,
                Contains,
            };

            struct GeneData
            {
                inline Locus l() const
                {
                    Base end   = std::numeric_limits<Base>::min();
                    Base start = std::numeric_limits<Base>::max();
                
                    for (const auto &i : seqs)
                    {
                        end   = std::max(end, i->l.end);
                        start = std::max(end, i->l.end);
                    }

                    return Locus(start, end);
                }

                // Calcualate the normalized abundance for the gene
                inline Concentration abund(Mixture m) const
                {
                    Concentration n = 0;
                
                    for (const auto &i : seqs)
                    {
                        n += ((*i).mixes.at(m) / (*i).length);
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
        
            // Return the detection limit at the gene level
            Sensitivity limitGene(const GeneHist &) const;

            // Number of non-overlapping bases in all exons
            Base exonBase() const;

            const std::vector<ExonData> & mergedExons()     const;
            const std::vector<ExonData> & sortedExons()     const;
            const std::vector<IntronData> & sortedIntrons() const;

            const GeneData   *findGene  (const GeneID &id)         const;
            const GeneData   *findGene  (const Locus &l, Matching) const;
            const ExonData   *findExon  (const Locus &l, Matching) const;
            const IntronData *findIntron(const Locus &l, Matching) const;

        private:

            struct TransRefImpl;

            std::shared_ptr<TransRefImpl> _impl;        
    };
}

#endif