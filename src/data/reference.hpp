#ifndef GI_REFERENCE_HPP
#define GI_REFERENCE_HPP

#include <memory>

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

    struct AnnotationData
    {
        AnnotationData(const SequinID &id, Locus l) : id(id), l(l) {}
        
        inline bool operator<(const AnnotationData &x)  const { return id < x.id;  }
        inline bool operator==(const AnnotationData &x) const { return id == x.id; }
        
        SequinID id;
        Locus l;
    };

    template <typename Data = SequinData, typename Stats = SequinStats, typename Annot = AnnotationData>
        class Reference
    {
        public:
            Reference()
            {
                _mixes  = std::shared_ptr<Mixtures>(new Mixtures);
                _annots = std::shared_ptr<Annotations>(new Annotations);            
            }

            // Add a sequin defined in a mixture file
            void add(const SequinID &id, Base length, Concentration c, Mixture m)
            {
                (*_mixes)[m].insert(MixtureData(id, length, c));
                _rawMIDs.insert(id);
            }

            // Return number of sequins in the mixture
            inline std::size_t countMixes() const { return _mixes->size(); }

            // Return all validated sequins
            inline const std::map<SequinID, Data> &data() const { return _data; }

            inline const Data *seq(const SequinID &id) const
            {
                return _data.count(id) ? &_data.at(id) : nullptr;
            }

            virtual void validate()
            {
//                std::vector<SequinID> mixIDs, antIDs;
//                
//                antIDs.resize((*_annots).size());
//                mixIDs.resize((*_mixes)[MixA].size());
//                
//                std::transform((*_mixes)[MixA].begin(), (*_mixes)[MixA].end(), mixIDs.begin(),
//                               [&](const MixtureData &m)
//                               {
//                                   return m.id;
//                               });
//                
//                std::transform((*_annots).begin(), (*_annots).end(), antIDs.begin(),
//                               [&](const AnnotationData &m)
//                               {
//                                   return m.id;
//                               });
//                
//                /*
//                 * Check for any sequin defined in mixture but not in annotation
//                 */
//                
//                std::vector<SequinID> diffs, inters;
//                
//                std::set_difference(mixIDs.begin(),
//                                    mixIDs.end(),
//                                    antIDs.begin(),
//                                    antIDs.end(),
//                                    std::back_inserter(diffs));
//                
//                std::set_intersection(mixIDs.begin(),
//                                      mixIDs.end(),
//                                      antIDs.begin(),
//                                      antIDs.end(),
//                                      std::back_inserter(inters));
//                
//                /*
//                 * Construct a set of validated sequins
//                 */
//                
//                std::for_each(mixIDs.begin(), mixIDs.end(), [&](const SequinID &id)
//                              {
//                                  auto d = SequinData();
//                                  
//                                  // The rest of the fields will be filled later...
//                                  d.id = id;
//                                  
//                                  // Add a new entry for the validated sequin
//                                  _data[id] = d;
//                              });
//                
//                /*
//                 * Now, we have a list of validated sequins. Use those sequins to combine information
//                 * from mixtures and annotations.
//                 */
//                
//                for (const auto i : (*_mixes))
//                {
//                    // Eg: MixA, MixB etc
//                    const auto mix = i.first;
//                    
//                    // For each of the mixture defined
//                    for (const auto j : i.second)
//                    {
//                        // Only if it's a validated sequin
//                        if (_data.count(j.id))
//                        {
//                            _data.at(j.id).length = j.length;
//                            _data.at(j.id).mixes[mix] = j.abund;
//                        }
//                    }
//                }
            }

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

            struct Mixtures : public std::map<Mixture, std::set<MixtureData>>
            {
                // Empty Implementation
            };

            struct Annotations : std::set<AnnotationData>
            {
                // Empty Implementation
            };

            std::set<SequinID> _rawMIDs;
            std::set<SequinID> _rawAIDs;

            // Validated sequins
            std::map<SequinID, Data> _data;

            // Statistics about the sequins, only valid after validate()
            Stats _stats;

            std::shared_ptr<Mixtures> _mixes;
            std::shared_ptr<Annotations> _annots;
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

//        inline ExonData operator+(const ExonData &l)  const
  //      {
    //        return ExonData(std::min(start, l.start), std::max(end, l.end));
      //  }

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

    class TransReference : public Reference<SequinID, SequinStats, ExonData>
    {
        public:
            inline void adds(const IsoformID &iID, const GeneID &gID, const Locus &l)
            {
                ExonData e;
            
                e.l   = l;
                e.iID = iID;
                e.gID = gID;
            
                _exonsByGenes[gID].push_back(e);
                _exonsByTrans[iID].push_back(e);
            
                _rawAIDs.insert(iID);
                _rawGIDs.insert(gID);
            }
        
            template <typename Data> std::map<SequinID, Data> merge(const std::set<SequinID> &mIDs, const std::set<SequinID> &aIDs)
            {
                assert(!mIDs.empty() && !aIDs.empty());
            
                /*
                 * Check for any sequin defined in mixture but not in annotation
                 */
            
                std::vector<SequinID> diffs, inters;
            
                std::set_difference(mIDs.begin(),
                                    mIDs.end(),
                                    aIDs.begin(),
                                    aIDs.end(),
                                    std::back_inserter(diffs));

                std::set_intersection(mIDs.begin(),
                                      mIDs.end(),
                                      aIDs.begin(),
                                      aIDs.end(),
                                      std::back_inserter(inters));
            
                std::map<SequinID, Data> merged;
            
                /*
                 * Construct a set of validated sequins
                 */
            
                std::for_each(inters.begin(), inters.end(), [&](const SequinID &id)
                {
                    auto d = SequinData();
                              
                    // The rest of the fields will be filled later...
                    d.id = id;
                              
                    // Add a new entry for the validated sequin
                    merged[id] = d;
                });
            
                /*
                 * Now, we have a list of validated sequins. Use those sequins to combine information
                 * from mixtures and annotations.
                 */
            
                for (const auto i : (*_mixes))
                {
                    // Eg: MixA, MixB etc
                    const auto mix = i.first;
                
                    // For each of the mixture defined
                    for (const auto j : i.second)
                    {
                        // Only if it's a validated sequin
                        if (_data.count(j.id))
                        {
                            merged.at(j.id).length = j.length;
                            merged.at(j.id).mixes[mix] = j.abund;
                        }
                    }
                }
            
                return merged;
            }
        
            inline void validate()
            {
                if (_exonsByGenes.empty())
                {
                    throw std::runtime_error("There is no synthetic chromosome in the annotation file. Anaquin is unable to proceed unless a valid annotation is given. Please check your file and try again.");
                }
            
                const auto r = merge<SequinData>(_rawMIDs, _rawAIDs);
            
                /*
                 * Filter out only those validated exons
                 */
            
                for (const auto &i : _exonsByTrans)
                {
                    if (r.count(i.first))
                    {
                        for (const auto &j : i.second)
                        {
                            _sortedExons.push_back(j);
                        }
                    }
                }
            
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
                
                _exonBase = countLocus(_mergedExons = Locus::merge<ExonData, ExonData>(_sortedExons));

            //            /*
            //             * Construct data-structure for each gene
            //             */
            //
            //            for (const auto &geneID : _exonsByGenes)
            //            {
            //                Feature g;
            //
            //                // The name of the gene is also the it's ID
            //                g.id = g.geneID = geneID;
            //
            //                g.l.end   = std::numeric_limits<Base>::min();
            //                g.l.start = std::numeric_limits<Base>::max();
            //
            //                /*
            //                 * Add all exons for this gene
            //                 */
            //
            //                for (const auto &f : fs)
            //                {
            //                    if (f.geneID == g.id)
            //                    {
            //                        g.l.end   = std::max(g.l.end, f.l.end);
            //                        g.l.start = std::min(g.l.start, f.l.start);
            //
            //                        if (f.type == Exon)
            //                        {
            //                            // TODO: ????
            //                            r_l_exons.push_back(RNALocus(f.geneID, f.l));
            //                        }
            //                    }
            //                }
            //
            //                assert(g.l.end   != std::numeric_limits<Base>::min());
            //                assert(g.l.start != std::numeric_limits<Base>::min());
            //                assert(g.l.end > g.l.start);
            //
            //                r_genes.push_back(g);
            //            }
            }
        
            inline Base exonBase() const { return _exonBase; }
        
            inline const std::vector<ExonData> & mergedExons()     const { return _mergedExons;   }
            inline const std::vector<ExonData> & sortedExons()     const { return _sortedExons;   }
            inline const std::vector<IntronData> & sortedIntrons() const { return _sortedIntrons; }

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

            inline const IntronData *findIntrons(const Locus &l) const
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
        
            std::vector<ExonData>   _mergedExons;
            std::vector<ExonData>   _sortedExons;
            std::vector<IntronData> _sortedIntrons;
        
            std::set<std::string>                      _rawGIDs;
            std::map<GeneID, std::vector<ExonData>>    _exonsByGenes;
            std::map<IsoformID, std::vector<ExonData>> _exonsByTrans;
    };
}

#endif