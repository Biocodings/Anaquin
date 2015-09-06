#ifndef GI_REFERENCE_HPP
#define GI_REFERENCE_HPP

#include <map>
#include "data/types.hpp"
#include "data/locus.hpp"

namespace Anaquin
{
    template <typename Iter> Base countLocus(const Iter &iter)
    {
        Base n = 0;
        
        for (const auto &i : iter)
        {
            n += static_cast<Locus>(i).length();
        }
        
        return n;
    }
    
    enum Mixture
    {
        MixA,
        MixB,
        MixF,
        MixG,
    };

    struct SequinStats
    {
    };

    struct SequinData
    {
        inline bool operator<(const SequinID &x)  const { return this->id < x;  }
        inline bool operator==(const SequinID &x) const { return this->id == x; }

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
        
            inline SequinHist hist() const
            {
                SequinHist h;
            
                for (const auto &i : _data)
                {
                    h[i.first] = 0;
                }

                return h;
            }

            virtual void validate() {};

        protected:

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

    struct TransSequinData : public SequinData
    {
        GeneID gID;
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
        
        inline Base length() const
        {
            Base n = 0;
            
            for (const auto &i : seqs)
            {
                n = std::max(static_cast<Base>(0), i->l.length());
            }

            return n;
        }
        
        inline Concentration abund(Mixture m = MixA) const
        {
            Concentration n = 0;
            
            for (const auto &i : seqs)
            {
                n += (*i).mixes.at(m);
            }

            return n;
        }
        
        // Each sequin comprises an isoform
        std::vector<SequinData *> seqs;
    };

    class TransRef : public Reference<TransSequinData, SequinStats>
    {
        public:

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

            // Return a histogram for all the validated genes
            inline GeneHist histForGene() const
            {
                GeneHist h;
            
                for (const auto &i : _data)
                {
                    h[i.second.gID] = 0;
                }

                return h;
            }

            // Add a new annoation reference
            inline void adds(const IsoformID &iID, const GeneID &gID, const Locus &l)
            {
                ExonData e;
            
                e.l   = l;
                e.iID = iID;
                e.gID = gID;
            
                _exonsByGenes[gID].push_back(e);
                _exonsByTrans[iID].push_back(e);
            
                _rawIIDs.insert(iID);
                _rawGIDs.insert(gID);
            }

            template <typename Data> void merge(const std::set<SequinID> &mIDs,
                                                const std::set<SequinID> &aIDs)
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
                    auto d = Data();

                    d.id = id;

                    // Add a new entry for the validated sequin
                    _data[id] = d;
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
                
                for (const auto &i : _exonsByTrans)
                {
                    if (_data.count(i.first))
                    {
                        _data[i.first].l = Locus::expand(i.second, [&](const ExonData &f)
                        {
                            return true;
                        });

                        _genes[_data[i.first].gID].seqs.push_back(&_data[i.first]);
                    }
                }
            }

            inline void validate()
            {
                if (_exonsByGenes.empty())
                {
                    throw std::runtime_error("There is no synthetic chromosome in the annotation file. Anaquin is unable to proceed unless a valid annotation is given. Please check your file and try again.");
                }

                merge<TransSequinData>(_rawMIDs, _rawIIDs);

                /*
                 * Filter out only those validated exons
                 */
                
                for (const auto &i : _exonsByTrans)
                {
                    if (_data.count(i.first))
                    {
                        for (const auto &j : i.second)
                        {
                            //_gIDs.insert(j.gID);
                            _sortedExons.push_back(j);
                        }
                    }
                }

                //assert(!_gIDs.empty());
        
                /*
                 * Generate a list of sorted exons
                 */
            
                assert(!_sortedExons.empty());
                std::sort(_sortedExons.begin(), _sortedExons.end(), [](const ExonData &x, const ExonData &y)
                {
                    return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
                });
            
                /*
                 * Generate a list of sorted introns, only possible once the exons are sorted.
                 */
            
                for (auto i = 0; i < _sortedExons.size(); i++)
                {
                    if (i && _sortedExons[i].gID == _sortedExons[i-1].gID)
                    {
                        IntronData d;
                    
                        d.gID = _sortedExons[i].gID;
                        d.iID = _sortedExons[i].iID;
                        d.l   = Locus(_sortedExons[i-1].l.end + 1, _sortedExons[i].l.start - 1);
                    
                        _sortedIntrons.push_back(d);
                    }
                }
            
                assert(!_sortedIntrons.empty());
                std::sort(_sortedIntrons.begin(), _sortedIntrons.end(), [](const IntronData &x, const IntronData &y)
                {
                    return (x.l.start < y.l.start) || (x.l.start == y.l.start && x.l.end < y.l.end);
                });
                
                // Count number of non-overlapping bases for all exons
                _exonBase = countLocus(_mergedExons = Locus::merge<ExonData, ExonData>(_sortedExons));
            }

            // Number of non-overlapping bases in all exons
            inline Base exonBase() const { return _exonBase; }
        
            inline const std::vector<ExonData> & mergedExons()     const { return _mergedExons;   }
            inline const std::vector<ExonData> & sortedExons()     const { return _sortedExons;   }
            inline const std::vector<IntronData> & sortedIntrons() const { return _sortedIntrons; }

            inline const GeneData *findGene(const GeneID &id) const
            {
                return _genes.count(id) ? &(_genes.at(id)) : nullptr;
            }

            inline const GeneData *findGene(const Locus &l) const
            {
                for (const auto &i : _genes)
                {
                    if (i.second.l().overlap(l))
                    {
                        return &i.second;
                    }
                }

                return nullptr;
            }

            inline const ExonData *findExon(const Locus &l) const
            {
                for (const auto &i : _sortedExons)
                {
                    if (i.l == l)
                    {
                        return &i;
                    }
                }
                
                return nullptr;
            }

            inline const IntronData *findIntron(const Locus &l) const
            {
                for (const auto &i : _sortedIntrons)
                {
                    if (i.l == l)
                    {
                        return &i;
                    }
                }
            
                return nullptr;
            }

        private:

            // Number of bases for all the reference exons
            Base _exonBase;

            //std::set<GeneID> _gIDs;
        
            std::map<GeneID, GeneData> _genes;
            std::vector<ExonData>   _mergedExons;
            std::vector<ExonData>   _sortedExons;
            std::vector<IntronData> _sortedIntrons;

            /*
             * Raw data - structure before validated
             */

            std::set<GeneID>                           _rawGIDs;
            std::set<IsoformID>                        _rawIIDs;
            std::map<GeneID, std::vector<ExonData>>    _exonsByGenes;
            std::map<IsoformID, std::vector<ExonData>> _exonsByTrans;
    };
}

#endif