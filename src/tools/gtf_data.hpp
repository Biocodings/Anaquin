#ifndef GTF_DATA_HPP
#define GTF_DATA_HPP

#include "data/hist.hpp"
#include "data/intervals.hpp"
#include "tools/gtf_data.hpp"
#include "parsers/parser_gtf.hpp"

namespace Anaquin
{
    struct ExonData_
    {
        operator const Locus &() const { return l; }
        
        inline bool operator<(const  ExonData_ &x) const { return l < x.l;  }
        inline bool operator==(const ExonData_ &x) const { return l == x.l; }
        inline bool operator!=(const Locus &l)    const { return !operator==(l); }
        inline bool operator==(const Locus &l)    const
        {
            return this->l.start == l.start && this->l.end == l.end;
        }
        
        inline void operator+=(const ExonData_ &x)
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
    
    typedef ExonData_ IntronData_;
    
    struct TransData_
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
    
    struct ChrData
    {
        // Transcripts to Data
        std::map<TransID, TransData_> t2d;
        
        // Genes to Data
        std::map<GeneID, GeneData> g2d;
        
        // Genes to exons
        std::map<TransID, std::set<ExonData_>> t2e;
        
        // Genes to introns
        std::map<TransID, std::set<IntronData_>> t2i;
        
        Counts exons = 0;
        Counts intrs = 0;
    };
    
    struct GTFData : public std::map<ChrID, ChrData>
    {
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
        
        inline Counts countIntr() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countIntr(cID);
            });
        }
        
        inline Counts countGene(const ChrID &cID) const
        {
            return at(cID).g2d.size();
        }
        
        inline Counts countTrans(const ChrID &cID) const
        {
            return at(cID).t2d.size();
        }
        
        inline Counts countExon(const ChrID &cID) const
        {
            return at(cID).exons;
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
        
        inline Counts countExonsGen() const
        {
            return countExon() - countExonSyn();
        }
        
        inline Counts countIntrGen() const
        {
            return countIntr() - countIntrSyn();
        }

        inline Intervals<> intervals(const ChrID &cID) const
        {
            Intervals<> r;
            
            for (const auto &i : at(cID).g2d)
            {
                r.add(Interval(i.first, i.second.l));
            }

            return r;
        }
        
        // Returns total length of all genes for a chromosome
        inline Base countLen(const ChrID &cID) const
        {
            // Assuming the genes are non-overlapping
            return intervals(cID).stats().length;
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
                    ExonData_ d;
                    
                    d.l   = x.l;
                    d.cID = x.cID;
                    d.gID = x.gID;
                    d.tID = x.tID;
                    
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
            for (const auto &j : i.second.t2e)
            {
                // Sorted exons
                auto sorted = std::vector<ExonData_>();
                
                // Convert to a sorted vector
                std::copy(j.second.begin(), j.second.end(), std::back_inserter(sorted));
                
                /*
                 * Generate a list of sorted introns, only possible once the exons are sorted.
                 */
                
                for (auto j = 1; j < sorted.size(); j++)
                {
                    const auto &x = sorted[j-1];
                    const auto &y = sorted[j];
                    
                    IntronData_ d;
                    
                    d.gID = x.gID;
                    d.tID = x.tID;
                    d.cID = x.cID;
                    d.l   = Locus(x.l.end+1, y.l.start-1);
                    
                    c2d[x.cID].intrs++;
                    i.second.t2i[d.tID].insert(d);
                }
            }
        }
        
        return c2d;
    }
}

#endif