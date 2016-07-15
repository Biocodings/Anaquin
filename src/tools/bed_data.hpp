#ifndef BED_DATA_HPP
#define BED_DATA_HPP

#include "data/hist.hpp"
#include "data/minters.hpp"
#include "data/standard.hpp"
#include "data/intervals.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    typedef ParserBed::Data Data;

    struct BedChrData
    {
        std::map<std::string, Data> g2d;
    };
    
    struct BedData : public std::map<ChrID, BedChrData>
    {
        inline Counts countGene(const ChrID &cID) const
        {
            return at(cID).g2d.size();
        }
        
        inline Counts countBase(const ChrID &cID) const
        {
            Base b = 0;
            
            for (const auto &i : at(cID).g2d)
            {
                b += i.second.l.length();
            }

            assert(b);
            return b;
        }

        inline Counts countBaseSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return Standard::isSynthetic(cID) ? countBase(cID) : 0;
            });
        }

        inline Counts countBase() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return countBase(cID);
            });
        }
        
        inline Counts countBaseGen() const
        {
            return countBase() - countBaseSyn();
        }
        
        inline Counts countGene() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return countGene(cID);
            });
        }
        
        inline Counts countGeneSyn() const
        {
            return ::Anaquin::count(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return Standard::isSynthetic(cID) ? countGene(cID) : 0;
            });
        }

        inline Counts countGeneGen() const
        {
            return countGene() - countGeneSyn();
        }

        inline std::map<ChrID, Hist> hist() const
        {
            std::map<ChrID, Hist> r;
            
            for (const auto &i : *this)
            {
                for (const auto &j :i.second.g2d)
                {
                    r[i.first][j.first];
                }
            }
            
            return r;
        }
        
        inline Intervals<> inters(const ChrID &cID) const
        {
            Intervals<> r;
            
            for (const auto &i : at(cID).g2d)
            {
                r.add(Interval(i.first, i.second.l));
            }
            
            r.build();
            return r;
        }

        // Intervals for the genes
        inline std::map<ChrID, Intervals<>> inters() const
        {
            std::map<ChrID, Intervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = inters(i.first);
            }
            
            return r;
        }
        
        // Genomic regions mapped by chromosomes
        inline std::map<ChrID, Intervals<>> intersGen() const
        {
            std::map<ChrID, Intervals<>> r;
            
            for (const auto &i : *this)
            {
                if (Standard::isGenomic(i.first))
                {
                    r[i.first] = inters(i.first);
                }
            }
            
            return r;
        }
        
        // Synthetic regions mapped by chromosomes
        inline ID2Intervals intersSyn() const
        {
            ID2Intervals r;
            
            for (const auto &i : *this)
            {
                if (Standard::isSynthetic(i.first))
                {
                    r[i.first] = inters(i.first);
                }
            }
            
            return r;
        }

        inline MergedIntervals<> minters(const ChrID &cID) const
        {
            MergedIntervals<> r;
            
            for (const auto &i : at(cID).g2d)
            {
                #define BIN_SIZE 10000
                
                if (false && i.second.l.length() >= BIN_SIZE)
                {
                    Locus l;
                    
                    l.start = 1;
                    l.end = BIN_SIZE;
                    
                    while (l.end < i.second.l.length())
                    {
                        r.add(MergedInterval(toString(i.first) + "_" + toString(l.start), l));
                        
                        l.start = l.end + 1;
                        l.end = l.start + BIN_SIZE;
                    }
                }
                else
                {
                    r.add(MergedInterval(i.first, i.second.l));
                }
            }
            
            r.build();
            
            return r;
        }
        
        inline std::map<ChrID, MergedIntervals<>> minters() const
        {
            std::map<ChrID, MergedIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = minters(i.first);
            }
            
            return r;
        }
        
        inline std::map<ChrID, MergedIntervals<>> mintersGen() const
        {
            std::map<ChrID, MergedIntervals<>> r;
            
            for (const auto &i : *this)
            {
                if (Standard::isGenomic(i.first))
                {
                    r[i.first] = minters(i.first);
                }
            }
            
            return r;
        }
        
        inline std::map<ChrID, MergedIntervals<>> mintersSyn() const
        {
            std::map<ChrID, MergedIntervals<>> r;
            
            for (const auto &i : *this)
            {
                if (Standard::isSynthetic(i.first))
                {
                    r[i.first] = minters(i.first);
                }
            }
            
            return r;
        }
    };
    
    inline BedData bedData(const Reader &r)
    {
        BedData c2d;
        
        ParserBed::parse(r, [&](const ParserBed::Data &x, const ParserProgress &)
        {
            if (Standard::isSynthetic(x.cID))
            {
                const auto bID = baseID(x.name);
                
                if (c2d[x.cID].g2d.count(bID))
                {
                    return;
                }

                c2d[x.cID].g2d[bID] = x;
            }
            else
            {
                c2d[x.cID].g2d[x.name] = x;
            }
        });

        return c2d;
    }    
}

#endif
