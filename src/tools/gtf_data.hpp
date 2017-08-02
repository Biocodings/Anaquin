#ifndef GTF_DATA_HPP
#define GTF_DATA_HPP

#include "data/hist.hpp"
#include "tools/tools.hpp"
#include "data/dinters.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "parsers/parser_gtf.hpp"

// Defined in main.cpp
extern void printWarning(const std::string &);

namespace Anaquin
{
    struct ExonData
    {
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

        inline bool isForward()  const { return str == Strand::Forward;  }
        inline bool isBackward() const { return str == Strand::Backward; }
        
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        // Eg: ENST00000456328.2
        TransID tID;
        
        Strand str;
        
        Locus l;
    };
    
    struct TransData
    {
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        // Eg: ENST00000456328.2
        TransID tID;
        
        Locus l;
    };
    
    struct GeneData
    {
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        Locus l;
    };

    typedef ExonData IntronData;
    
    struct ChrData
    {
        // Transcripts to Data
        std::map<TransID, TransData> t2d;
        
        // Genes to Data
        std::map<GeneID, GeneData> g2d;
        
        // Transcripts to non-unique exons
        std::map<TransID, std::set<ExonData>> t2e;
        
        // Transcripts to unique exons
        std::map<TransID, std::set<ExonData>> t2ue;

        // Transcripts to unique introns
        std::map<TransID, std::set<IntronData>> t2ui;
        
        // Transcripts to genes
        std::map<TransID, GeneID> t2g;
        
        // Unique exons (forwards + backwards)
        Counts uexons = 0;

        // Unique introns (forwards + backwards)
        Counts uintrs = 0;
        
        std::set<GeneID> gIDs;
    };
    
    struct GTFData : public std::map<ChrID, ChrData>
    {
        // Genes for a chromosome
        inline std::set<GeneID> genes(const ChrID &x) const
        {
            std::set<GeneID> r;
            
            for (const auto &i : at(x).g2d)
            {
                r.insert(i.first);
            }
            
            return r;
        }
        
        // Genes for all chromosome
        inline std::map<ChrID, std::set<GeneID>> genes() const
        {
            std::map<ChrID, std::set<GeneID>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = genes(i.first);
            }
            
            return r;
        }
        
        inline Counts nGene() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return nGene(cID);
            });
        }
        
        inline Counts countTrans() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countTrans(cID);
            });
        }
        
        inline Counts countUExon() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countUExon(cID);
            });
        }

        inline Counts countUIntr() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countUIntr(cID);
            });
        }
        
        inline Counts nGene(const ChrID &cID) const
        {
            return at(cID).gIDs.size();
        }
        
        inline Counts countTrans(const ChrID &cID) const
        {
            return at(cID).t2d.size();
        }
        
        inline Counts countUExon(const ChrID &cID) const
        {
            return at(cID).uexons;
        }
        
        inline Counts countUIntr(const ChrID &cID) const
        {
            return at(cID).uintrs;
        }

        inline Counts nGeneSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return isChrIS(cID) ? nGene(cID) : 0;
            });
        }
        
        inline Counts countTransSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return isChrIS(cID) ? countTrans(cID) : 0;
            });
        }
        
        inline Counts countUExonSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return isChrIS(cID) ? countUExon(cID) : 0;
            });
        }
        
        inline Counts countUIntrSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return isChrIS(cID) ? countUIntr(cID) : 0;
            });
        }

        inline Counts nGeneGen() const
        {
            return nGene() - nGeneSyn();
        }
        
        inline Counts countTransGen() const
        {
            return countTrans() - countTransSyn();
        }
        
        inline Counts countUExonGen() const
        {
            return countUExon() - countUExonSyn();
        }

        inline Counts countUIntrGen() const
        {
            return countUIntr() - countUIntrSyn();
        }
        
        inline DIntervals<> gIntervals(const ChrID &cID) const
        {
            DIntervals<> r;
            
            for (const auto &i : at(cID).g2d)
            {
                r.add(DInter(i.first, i.second.l));
            }
            
            r.build();
            return r;
        }
        
        inline std::map<ChrID, DIntervals<>> gIntervals() const
        {
            std::map<ChrID, DIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = gIntervals(i.first);
            }

            return r;
        }
        
        // Intervals for merged exons
        inline Chr2MInters meInters(Strand str) const
        {
            Chr2MInters r;
            
            for (const auto &i : *this)
            {
                r[i.first] = meInters(i.first, str);
            }
            
            return r;
        }
        
        // Intervals for merged exons (only possible at the gene level)
        inline MergedIntervals<> meInters(const ChrID &cID, Strand str) const
        {
            MergedIntervals<> r;

            // This is needed to merge exons over all transcripts
            std::map<GeneID, MergedInterval> merged;
            
            // For each transcript...
            for (const auto &i : at(cID).t2ue)
            {
                const auto &tID = i.first;
                const auto &gID = at(cID).t2g.at(tID);

                if (!merged.count(gID))
                {
                    // Merging unique exons (doesn't matter how long this is)
                    merged[gID] = MergedInterval(gID, Locus(1, std::numeric_limits<Base>::max()));
                }
                
                // For each exon in the transcript...
                for (const auto &j : i.second)
                {
                    if (str != Strand::Either && (j.str != str))
                    {
                        break;
                    }

                    // Merge all the overlapping exons
                    merged.at(gID).map(j.l);
                }
            }
            
            // For each gene in the chromosome...
            for (const auto &i : merged)
            {
                const auto &gID = i.first;

                // For each merged exon in the gene...
                for (const auto &j : i.second._data)
                {
                    const auto &l = j.second;
                    r.add(MergedInterval(gID + "-" + toString(l.start) + "-" + toString(l.end), l, gID, gID));
                }
            }

            r.build();
            return r;
        }

        // Intervals for unique exons
        inline MergedIntervals<> ueInters(const ChrID &cID) const
        {
            MergedIntervals<> r;
            
            for (const auto &i : at(cID).t2ue)
            {
                const auto &tID = i.first;
                const auto &gID = at(cID).t2g.at(tID);
                
                for (const auto &j : i.second)
                {
                    r.add(MergedInterval(i.first + "-" + toString(j.l.start) + "-" + toString(j.l.end),
                                         j.l, gID, tID));
                }
            }

            r.build();
            return r;
        }

        // Intervals for unique exons
        inline Chr2MInters ueInters() const
        {
            Chr2MInters r;
            
            for (const auto &i : *this)
            {
                r[i.first] = ueInters(i.first);
            }
            
            return r;
        }

        /*
         * Returns non-overlapping intervals for a chromosome. Each interval represents a set of merged
         * overlapping unique exons.
         */
        
        inline MergedIntervals<> mergedExons(const ChrID &cID) const
        {
            MergedIntervals<> r;
            
            for (const auto &i : at(cID).t2ue)
            {
                const auto &tID = i.first;
                const auto &gID = at(cID).t2g.at(tID);

                for (auto &j : i.second)
                {
                    r.merge(MergedInterval(i.first + "-" + toString(j.l.start) + "-" + toString(j.l.end),
                                           j.l,
                                           gID,
                                           tID));
                }
            }

            r.build();
            return r;
        }
        
        /*
         * Returns non-overlapping exon intervals for all chromosomes.
         */

        inline Chr2MInters mergedExons() const
        {
            Chr2MInters r;
            
            for (const auto &i : *this)
            {
                r[i.first] = mergedExons(i.first);
            }
            
            return r;
        }
        
        // Intervals for unique introns
        inline MergedIntervals<> uiInters(const ChrID &cID) const
        {
            MergedIntervals<> r;
            
            for (const auto &i : at(cID).t2ui)
            {
                const auto &tID = i.first;
                const auto &gID = at(cID).t2g.at(tID);

                for (const auto &j : i.second)
                {
                    r.add(MergedInterval(i.first + "-" + toString(j.l.start) + "-" + toString(j.l.end),
                                         j.l, gID, tID));
                }
            }

            // Eg: chrM doesn't have any intron...
            if (r.size())
            {
                r.build();
            }
            
            return r;
        }
        
        // Intervals for unique introns
        inline Chr2MInters uiInters() const
        {
            Chr2MInters r;
            
            for (const auto &i : *this)
            {
                r[i.first] = uiInters(i.first);
            }
            
            return r;
        }
        
        // Returns total length of all genes for a chromosome
        inline Base countLen(const ChrID &cID) const
        {
            // Assuming the genes are non-overlapping
            return gIntervals(cID).stats().length;
        }
        
        inline Base countLenSyn() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return isChrIS(cID) ? countLen(cID) : 0;
            });
        }
        
        inline Base countLenGen() const
        {
            return countMap(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return !isChrIS(cID) ? countLen(cID) : 0;
            });
        }
    };

    inline GTFData gtfData(const Reader &r)
    {
        GTFData c2d;
        
        // Used for unique exons
        std::map<std::string, Locus> m_exons;
        
        // Used for unique introns
        std::map<std::string, Locus> m_intrs;

        ExonData ed;
        GeneData gd;
        TransData td;
        IntronData id;
        
        ParserGTF::parse(r, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &)
        {
            switch (x.type)
            {
                case RNAFeature::Transcript:
                {
                    td.l   = x.l;
                    td.cID = x.cID;
                    td.gID = x.gID;
                    td.tID = x.tID;

                    c2d[x.cID].gIDs.insert(td.gID);
                    c2d[x.cID].t2g[td.tID] = td.gID;
                    c2d[x.cID].t2d[td.tID] = td;
                    break;
                }

                case RNAFeature::Gene:
                {
                    gd.l   = x.l;
                    gd.cID = x.cID;
                    gd.gID = x.gID;
                    
                    c2d[x.cID].g2d[gd.gID] = gd;
                    break;
                }

                case RNAFeature::Exon:
                {
                    ed.l   = x.l;
                    ed.str = x.str;
                    ed.cID = x.cID;
                    ed.gID = x.gID;
                    ed.tID = x.tID;

                    c2d[x.cID].t2e[ed.tID].insert(ed);

                    // The key represents a unique exon
                    const auto key = (boost::format("%1%_%2%_%3%-%4%") % x.cID
                                                                       % x.str
                                                                       % x.l.start
                                                                       % x.l.end).str();

                    // Make sure it's unique due to alternative splicing
                    if (!m_exons.count(key))
                    {
                        m_exons[key] = x.l;
                        c2d[x.cID].uexons++;
                        c2d[x.cID].t2ue[ed.tID].insert(ed);
                    }
                    
                    break;
                }

                // Eg: CDS
                default: { break; }
            }
        });
        
        /*
         * The information we have is sufficient for exons, transcripts and genes. We just
         * need to compute introns.
         */
        
        // For each chromosome...
        for (auto &i : c2d)
        {
            // For each transcript...
            for (const auto &j : i.second.t2e)
            {
                // Sorted exons
                auto sorted = std::vector<ExonData>();
                
                // Sort the exons
                std::copy(j.second.begin(), j.second.end(), std::back_inserter(sorted));
                
                /*
                 * -------------------------------------- Cufflinks Bug --------------------------------------
                 *
                 * It's possible for Cufflink guided assembly to give a transcript in both forward and backward
                 * strand. An example is in "cufflink_bug.png" in the source distribution. Clearly invalid,
                 * and thus we'll ignore the transcript.
                 */
                
                bool shouldSkip = false;

                for (auto k = 1; k < sorted.size() && !shouldSkip; k++)
                {
                    const auto &x = sorted[k-1];
                    const auto &y = sorted[k];

                    if (x.str != y.str)
                    {
                        printWarning(j.first + " gives transcription in both forward and backward strand. Ignored.");
                        shouldSkip = true;
                    }
                }

                /*
                 * Generating introns, only possible once the exons are sorted.
                 */
                
                for (auto j = 1; j < sorted.size() && !shouldSkip; j++)
                {
                    const auto &x = sorted[j-1];
                    const auto &y = sorted[j];

                    id.cID = x.cID;
                    id.gID = x.gID;
                    id.tID = x.tID;
                    id.str = x.str;

                    // Intron spans between exons
                    id.l = Locus(x.l.end+1, y.l.start-1);

                    #define MIN_INTRON_LEN 4
                    
                    const auto key = (boost::format("%1%_%2%-%3%") % x.cID
                                                                   % id.l.start
                                                                   % id.l.end).str();
                    
                    if (isChrIS(i.first) || id.l.length() >= MIN_INTRON_LEN)
                    {
                        if (!m_intrs.count(key))
                        {
                            c2d[x.cID].uintrs++;
                            m_intrs[key] = id.l;
                            c2d[x.cID].t2ui[id.tID].insert(id);
                        }
                    }
                }
            }
        }

        return c2d;
    }
}

#endif
