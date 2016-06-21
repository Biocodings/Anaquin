#ifndef GTF_DATA_HPP
#define GTF_DATA_HPP

#include <fstream>
#include "data/hist.hpp"
#include "data/intervals.hpp"
#include "parsers/parser_gtf.hpp"

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
        
        // Eg: chr1
        ChrID cID;
        
        // Eg: ENSG00000223972.5
        GeneID gID;
        
        // Eg: ENST00000456328.2
        TransID tID;
        
        Locus l;
    };
    
    typedef ExonData IntronData;
    
    struct ChrData
    {
        // Transcripts to Data
        std::map<TransID, TransData_> t2d;
        
        // Genes to Data
        std::map<GeneID, GeneData> g2d;
        
        // Genes to non-unique exons
        std::map<TransID, std::set<ExonData>> t2e;
        
        // Genes to non-unique introns
        std::map<TransID, std::set<IntronData>> t2i;

        // Genes to unique exons
        std::map<TransID, std::set<ExonData>> t2ue;
        
        // Genes to unique introns
        std::map<TransID, std::set<IntronData>> t2ui;
        
        // Transcripts to genes
        std::map<TransID, GeneID> t2g;
        
        // Non-unique exons
        Counts exons = 0;
        
        // Unique exons
        Counts uexons = 0;

        // Non-unique introns
        Counts intrs = 0;
        
        // Unique introns
        Counts uintrs = 0;
        
        std::set<GeneID> gIDs;
    };
    
    struct GTFData : public std::map<ChrID, ChrData>
    {
        // Position for exons (used for determining unique exons)
        std::map<ChrID, std::map<Locus, Counts>> el;
        
        // Position for introns (usd for determining unique introns)
        std::map<ChrID, std::map<Locus, Counts>> il;

        inline Hist histGene(const ChrID &cID) const
        {
            return createHist(at(cID).g2d);
        }

        inline std::map<ChrID, Hist> histGene() const
        {
            std::map<ChrID, Hist> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = histGene(i.first);
            }

            return r;
        }
        
        inline Counts countGene() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countGene(cID);
            });
        }
        
        inline Counts countTrans() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countTrans(cID);
            });
        }
        
        inline Counts countExon() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countExon(cID);
            });
        }

        inline Counts countUExon() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countUExon(cID);
            });
        }

        inline Counts countUIntr() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countUIntr(cID);
            });
        }
        
        inline Counts countIntr() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countIntr(cID);
            });
        }
        
        inline Counts countGene(const ChrID &cID) const
        {
            return at(cID).gIDs.size();
        }
        
        inline Counts countTrans(const ChrID &cID) const
        {
            return at(cID).t2d.size();
        }
        
        inline Counts countExon(const ChrID &cID) const
        {
            return at(cID).exons;
        }

        inline Counts countUExon(const ChrID &cID) const
        {
            return at(cID).uexons;
        }
        
        inline Counts countUIntr(const ChrID &cID) const
        {
            return at(cID).uintrs;
        }

        inline Counts countIntr(const ChrID &cID) const
        {
            return at(cID).intrs;
        }
        
        inline Counts countGeneSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countGene(cID) : 0;
            });
        }
        
        inline Counts countTransSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countTrans(cID) : 0;
            });
        }
        
        inline Counts countExonSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countExon(cID) : 0;
            });
        }
        
        inline Counts countUExonSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countUExon(cID) : 0;
            });
        }
        
        inline Counts countUIntrSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countUIntr(cID) : 0;
            });
        }

        inline Counts countIntrSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countIntr(cID) : 0;
            });
        }
        
        inline Counts countGeneGen() const
        {
            return countGene() - countGeneSyn();
        }
        
        inline Counts countTransGen() const
        {
            return countTrans() - countTransSyn();
        }
        
        inline Counts countExonGen() const
        {
            return countExon() - countExonSyn();
        }
        
        inline Counts countUExonGen() const
        {
            return countUExon() - countUExonSyn();
        }

        inline Counts countUIntrGen() const
        {
            return countUIntr() - countUIntrSyn();
        }
        
        inline Counts countIntrGen() const
        {
            return countIntr() - countIntrSyn();
        }

        inline Intervals<> gIntervals(const ChrID &cID) const
        {
            Intervals<> r;
            
            for (const auto &i : at(cID).g2d)
            {
                r.add(Interval(i.first, i.second.l));
            }
            
            r.build();
            return r;
        }
        
        inline std::map<ChrID, Intervals<>> gIntervals() const
        {
            std::map<ChrID, Intervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = gIntervals(i.first);
            }

            return r;
        }
        
        // Intervals for non-unique exons
        inline Intervals<> eIntervals(const ChrID &cID) const
        {
            Intervals<> r;
            
            for (const auto &i : at(cID).t2e)
            {
                for (const auto &j : i.second)
                {
                    // Each interval is an non-unique exon
                    r.add(Interval(i.first + "-" + toString(j.l.start) + "-" + toString(j.l.end), j.l));
                }
            }

            r.build();
            return r;
        }

        // Intervals for non-unique exons        
        inline std::map<ChrID, Intervals<>> eIntervals() const
        {
            std::map<ChrID, Intervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = eIntervals(i.first);
            }
            
            return r;
        }

        // Intervals for non-unique introns
        inline Intervals<> iIntervals(const ChrID &cID) const
        {
            Intervals<> r;
            
            for (const auto &i : at(cID).t2i)
            {
                for (const auto &j : i.second)
                {
                    // Each interval is an non-unique intron
                    r.add(Interval(i.first + "-" + toString(j.l.start) + "-" + toString(j.l.end), j.l));
                }
            }
            
            r.build();
            return r;
        }
        
        // Intervals for non-unique introns
        inline std::map<ChrID, Intervals<>> iIntervals() const
        {
            std::map<ChrID, Intervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = iIntervals(i.first);
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
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countLen(cID) : 0;
            });
        }
        
        inline Base countLenGen() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return !Standard::isSynthetic(cID) ? countLen(cID) : 0;
            });
        }
    };

    inline GTFData gtfData(const Reader &r)
    {
        GTFData c2d;
        
        ParserGTF::parse(r, [&](const ParserGTF::Data &x, const std::string &, const ParserProgress &)
        {
            switch (x.type)
            {
                case Transcript:
                {
                    TransData_ d;
                    
                    d.l   = x.l;
                    d.cID = x.cID;
                    d.gID = x.gID;
                    d.tID = x.tID;

                    c2d[x.cID].gIDs.insert(x.gID);
                    c2d[x.cID].t2g[d.tID] = d.gID;
                    c2d[x.cID].t2d[d.tID] = d;
                    break;
                }
                    
                case Gene:
                {
                    GeneData d;
                    
                    d.l   = x.l;
                    d.cID = x.cID;
                    d.gID = x.gID;
                    
                    c2d[x.cID].g2d[d.gID] = d;
                    break;
                }
                    
                case Exon:
                {
                    ExonData d;
                    
                    d.l   = x.l;
                    d.cID = x.cID;
                    d.gID = x.gID;
                    d.tID = x.tID;
                    
                    c2d.el[x.cID][d.l]++;
                    
                    c2d[x.cID].exons++;
                    c2d[x.cID].t2e[d.tID].insert(d);
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
        
        for (auto &i : c2d)
        {
            // For each transcript...
            for (const auto &j : i.second.t2e)
            {
                // Sorted exons
                auto sorted = std::vector<ExonData>();
                
                // Convert to a sorted vector
                std::copy(j.second.begin(), j.second.end(), std::back_inserter(sorted));
                
                /*
                 * Generate a list of sorted introns, only possible once the exons are sorted.
                 */

                for (auto j = 1; j < sorted.size(); j++)
                {
                    const auto &x = sorted[j-1];
                    const auto &y = sorted[j];
                    
                    IntronData d;
                    
                    d.gID = x.gID;
                    d.tID = x.tID;
                    d.cID = x.cID;
                    
                    // Intron simply spans between exons
                    d.l = Locus(x.l.end+1, y.l.start-1);

                    #define MIN_INTRON_LEN 20
                    
                    if (Standard::isSynthetic(i.first) || d.l.length() > MIN_INTRON_LEN)
                    {
                        c2d.il[x.cID][d.l]++;
                        c2d[x.cID].intrs++;
                        i.second.t2i[d.tID].insert(d);
                    }
                }
            }
        }
        
        /*
         * Calculating number of unique exons
         */

        for (auto &i : c2d.el)
        {
            c2d.at(i.first).uexons = i.second.size();
        }

        /*
         * Calculating number of unique introns
         */
        
        for (auto &i : c2d.il)
        {
//            if (i.first == "chr21")
//            {
//                std::ofstream out("/Users/tedwong/Sources/QA/data/RnaQuin/myChr21.txt");
//                
//                for (auto &j : i.second)
//                {//chr21   5012686 5014385
//                    out << "chr21\t" << j.first.start-2 << "\t" << j.first.end << std::endl;
//                }
//                
//                out.close();
//            }
            
            c2d.at(i.first).uintrs = i.second.size();
        }

        return c2d;
    }
}

#endif
