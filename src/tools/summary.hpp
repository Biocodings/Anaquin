#ifndef SUMMARY_HPP
#define SUMMARY_HPP

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
    
    struct Data : public std::map<ChrID, ChrData>
    {
        inline Counts countGenes() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countGenes(cID);
            });
        }
        
        inline Counts countTrans() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countTrans(cID);
            });
        }
        
        inline Counts countExons() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countExons(cID);
            });
        }
        
        inline Counts countIntrs() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return countIntrs(cID);
            });
        }
        
        inline Counts countGenes(const ChrID &cID) const
        {
            return at(cID).g2d.size();
        }
        
        inline Counts countTrans(const ChrID &cID) const
        {
            return at(cID).t2d.size();
        }
        
        inline Counts countExons(const ChrID &cID) const
        {
            return at(cID).exons;
        }
        
        inline Counts countIntrs(const ChrID &cID) const
        {
            return at(cID).intrs;
        }
        
        inline Counts countGenesSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countGenes(cID) : 0;
            });
        }
        
        inline Counts countTransSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countTrans(cID) : 0;
            });
        }
        
        inline Counts countExonsSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countExons(cID) : 0;
            });
        }
        
        inline Counts countIntrsSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const ChrData &x)
            {
                return Standard::isSynthetic(cID) ? countIntrs(cID) : 0;
            });
        }
        
        inline Counts countGenesGen() const
        {
            return countGenes() - countGenesSyn();
        }
        
        inline Counts countTransGen() const
        {
            return countTrans() - countTransSyn();
        }
        
        inline Counts countExonsGen() const
        {
            return countExons() - countExonsSyn();
        }
        
        inline Counts countIntrsGen() const
        {
            return countIntrs() - countIntrsSyn();
        }
    };
    
    Data c2d;
    
    static Data summaryGTF(const Reader &r)
    {
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
                    
                    c2d[x.gID].exons++;
                    c2d[x.gID].t2e[d.tID].insert(d);
                    break;
                }
                    
                    // Eg: CDS
                default: { break; }
            }
        });
        
        /*
         * The information we have is sufficient for answering exons, transcripts and genes. We just
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
                    
                    c2d[x.gID].intrs++;
                    i.second.t2i[d.tID].insert(d);
                }
            }
        }
        
        return c2d;
    }
}

#endif