#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include <set>
#include <map>
#include "stats/limit.hpp"
#include "data/variant.hpp"
#include "data/intervals.hpp"

namespace Anaquin
{
    typedef std::map<SequinID, Counts> EndoHist;
    typedef std::map<SequinID, Counts> SequinHist;

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
     * Generic template for a sequin. The definitions are appropriate for all sequins. Specalized
     * sequins are expected to derive from this class.
     */
    
    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

        // Return the abundance for this sequin specified by the mixture        
        inline Concent abund(Mixture m, bool norm = false) const
        {
            return mixes.at(m) / (norm ? l.length() : 1);
        }

        SequinID id;

        // Spiked-in concentration (not available if no mixture provided)
        std::map<Mixture, Concent> mixes;

        Locus l;
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

            template <typename T> SequinHist hist(const T &data) const
            {
                SequinHist h;
            
                for (const auto &i : data)
                {
                    h[i.first] = 0;
                }

                return h;
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
        
            // Calculate the detection limits
            inline Limit limit(const SequinHist &h) const
            {
                return limit(h, [&](const SequinID &id)
                {
                    return this->match(id);
                });
            }

        protected:

            virtual void validate() = 0;

            /*
             * Provides a common mechanism to calculate limit of detection given a histogram or a
             * distribution.
             */

            template <typename F> Limit limit
                                (const std::map<std::string, Counts> &h, F f, Mixture m = Mix_1) const
            {
                Limit s;
            
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
            void abund(const JoinID &, Concent &, Concent &, Concent &, Concent &, Mixture) const;

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
            
                // The position of the break-point
                Base l1, l2;
            
                // Orientation for each of the segment
                Strand s1, s2;
            };

            /*
             * This class represents the concentration between normal and fusion gene
             */

            struct SpliceChimeric
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
        
            // Add splicing for an intron in the normal gene
            void addSplice(const SequinID &, const Locus &);

            // Add fusion or normal genes comprised of the standards
            void addStand(const SequinID &, const Locus &);

            /*
             * Query operations
             */

            // Histogram for each known fusion
            SequinHist fusionHist() const;

            Counts countFusion() const;
            Counts countSplice() const;

            const KnownFusion *findFusion(const SequinID &) const;
            const KnownFusion *findFusion(Base x, Base y, Strand o1, Strand o2, double fuzzy) const;

            const SequinData *findSplice(const Locus &) const;
            const SpliceChimeric *findSpliceChim(const SequinID &) const;

            // Convert the normal gene to it's equivalent fusion (eg: NG1_1_P1 to FG1_1_P1)
            SequinID normalToFusion(const SequinID &) const;

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

            VarRef();

            /*
             * Modifier functions
             */
        
            // Adds a known variant
            void addVar(const Variant &);

            // Adds a sequin in the standards
            void addStand(const SequinID &, const Locus &);

            // Adds a reference interval (eg: chr21)
            void addRInterval(const ChromoID &, const Interval &);

            /*
             * Query functions
             */

            ChromoID endoID() const;

            bool isEndoID(const ChromoID &cID) const { return cID == endoID(); }

            const Intervals<> endoInters() const;
        
            // Returns number of known variants
            Counts countVars() const;

            // Returns number of known SNPs
            Counts countSNPs() const;
        
            // Returns number of known indels
            Counts countIndels() const;

            // Returns number of sequins
            Counts countSeqs() const;

            // Returns number of reference intervals
            Counts countInters() const;

            EndoHist endoHist() const;

            const Variant *findVar(const SequinID &) const;
            const Variant *findVar(const Locus &, MatchRule = Exact) const;

            Interval *findEndo(const Locus &) const;
            Interval *findEndo(const ChromoID &cID, const Locus &l) const
            {
                return isEndoID(cID) ? findEndo(l) : nullptr;
            }

            // Returns the expected fold-change
            Fold fold(const SequinID &) const;
        
            // Returns the expected allele frequency
            Proportion alleleFreq(const SequinID &) const;

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
    
    struct TransData : public SequinData
    {
        GeneID gID;
    };

    class TransRef : public Reference<TransData, DefaultStats>
    {
        public:

            struct ExonInterval : public Interval
            {
                ExonInterval(const GeneID    &gID,
                             const IsoformID &iID,
                             const IntervalID &id,
                             const Locus &l) : Interval(id, l), gID(gID), iID(iID) {}

                const GeneID    gID;
                const IsoformID iID;
            };
        
            typedef ExonInterval IntronInterval;
        
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
                inline Concent abund(Mixture m) const
                {
                    Concent n = 0;
                
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
                ExonData(const ChromoID  &cID,
                         const IsoformID &iID,
                         const GeneID    &gID,
                         const Locus &l)
                                : cID(cID), iID(iID), gID(gID), l(l) {}

                operator const Locus &() const { return l; }

                inline bool operator<(const  ExonData &x) const { return l < x.l;  }
                inline bool operator==(const ExonData &x) const { return l == x.l; }
                inline bool operator!=(const Locus &l)    const { return !operator==(l); }
                inline bool operator==(const Locus &l)    const
                {
                    return this->l.start == l.start && this->l.end == l.end;
                }
            
                inline void operator+=(const ExonData &x)
                {
                    l.start = std::min(l.start, x.l.start);
                    l.end   = std::max(l.end, x.l.end);
                }
            
                Locus     l;
                GeneID    gID;
                ChromoID  cID;
                IsoformID iID;
            };

            struct IntronData
            {
                operator const Locus &() const { return l; }
            
                inline bool operator<(const ExonData &x)  const { return l < x.l;  }
                inline bool operator==(const ExonData &x) const { return l == x.l; }
            
                Locus     l;
                GeneID    gID;
                ChromoID  cID;
                IsoformID iID;
            };

            TransRef();

            SequinHist geneHist(const ChromoID &) const;

            // Intervals for reference exons
            Intervals<ExonInterval> exonInters(const ChromoID &) const;
        
            // Intervals for reference introns
            Intervals<IntronInterval> intronInters(const ChromoID &) const;
        
            void addGene(const ChromoID &, const GeneID    &, const Locus &);
            void addExon(const ChromoID &, const IsoformID &, const GeneID &, const Locus &);

            /*
             * Accessor functions
             */

            Limit limitExon(const SequinHist &) const;
            Limit limitGene(const SequinHist &) const;
            Limit limitIsof(const SequinHist &) const;

            // Number of non-overlapping bases in all exons
            Base exonBase(const ChromoID &) const;

            // Number of merged exons
            Counts countMerged(const ChromoID &) const;
        
            // Number of unmerged exons
            Counts countExons(const ChromoID &) const;
        
            // Return number of introns
            Counts countIntrons(const ChromoID &) const;

            // List of chromosomes in the reference
            std::set<ChromoID> chromoIDs() const;
        
            inline bool isOnlyChrT() const { return chromoIDs().size() == 1; }

            std::vector<GeneID> geneIDs(const ChromoID &) const;

            const GeneData   *findGene  (const ChromoID &, const GeneID &)           const;
            const GeneData   *findGene  (const ChromoID &, const Locus &, MatchRule) const;
            const ExonData   *findExon  (const ChromoID &, const Locus &, MatchRule) const;
            const IntronData *findIntron(const ChromoID &, const Locus &, MatchRule) const;

        protected:
        
            void validate() override;

            void merge(const std::set<SequinID> &, const std::set<SequinID> &);
        
        private:

            struct TransRefImpl;

            std::shared_ptr<TransRefImpl> _impl;        
    };
}

#endif