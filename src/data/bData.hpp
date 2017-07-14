#ifndef B_DATA_HPP
#define B_DATA_HPP

#include "data/data.hpp"
#include "tools/tools.hpp"
#include "tools/errors.hpp"
#include "data/minters.hpp"
#include "data/dinters.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    typedef ParserBed::Data Data;

    struct BedChrData
    {
        std::map<SequinID, Data> r2d;
    };
    
    struct BedData : public std::map<ChrID, BedChrData>
    {
        inline std::set<SequinID> seqs() const
        {
            std::set<SequinID> x;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.r2d)
                {
                    x.insert(j.first);
                }
            }
            
            return x;
        }

        inline Counts nGene(const ChrID &cID) const
        {
            return at(cID).r2d.size();
        }

        inline Counts lengthReg(const SequinID &rID) const
        {
            for (const auto &i : *this)
            {
                if (i.second.r2d.count(rID))
                {
                    return i.second.r2d.at(rID).l.length();
                }
            }
            
            A_THROW("Unable to find " + rID);
        }

        inline Counts countBase(const ChrID &cID) const
        {
            Base b = 0;
            
            for (const auto &i : at(cID).r2d)
            {
                b += i.second.l.length();
            }

            A_ASSERT(b);
            return b;
        }

        inline Counts length() const
        {
            return countMap(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return countBase(cID);
            });
        }
        
        inline Counts nGene() const
        {
            return countMap(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return nGene(cID);
            });
        }

        inline Counts count() const
        {
            return countMap(*this, [&](const ChrID &, const BedChrData &x)
            {
                return x.r2d.size();
            });
        }
        
        template <typename F> Counts nGeneSyn(F f) const
        {
            return countMap(*this, [&](const ChrID &cID, const BedChrData &x)
            {
                return f(cID) ? nGene(cID) : 0;
            });
        }

        template <typename F> Counts nGeneGen(F f) const
        {
            return nGene() - nGeneSyn(f);
        }

        inline DIntervals<> inters(const ChrID &cID) const
        {
            DIntervals<> r;
            
            for (const auto &i : at(cID).r2d)
            {
                r.add(DInter(i.first, i.second.l));
            }
            
            r.build();
            return r;
        }

        inline std::map<ChrID, DIntervals<>> inters() const
        {
            std::map<ChrID, DIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = inters(i.first);
            }
            
            return r;
        }
        
        inline MergedIntervals<> minters(const ChrID &cID) const
        {
            MergedIntervals<> r;
            
            for (const auto &i : at(cID).r2d)
            {
                #define BIN_SIZE 10000
                
                if (false && i.second.l.length() >= BIN_SIZE)
                {
                    Locus l;
                    
                    l.start = 1;
                    l.end = BIN_SIZE;
                    
                    while (l.end < i.second.l.length())
                    {
                        r.add(MergedInterval(i.first + "_" + toString(l.start), l));
                        
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
    };
    
    struct RegionOptions
    {
        Base trim = 0;
    };
    
    template <typename F> BedData readRegions(const Reader &r, F f, RegionOptions o = RegionOptions())
    {
        BedData c2d;
        
        ParserBed::parse(r, [&](ParserBed::Data &x, const ParserProgress &p)
        {
            if (x.l.length() < 2 * o.trim)
            {
                throw std::runtime_error(x.name + " is too narrow. You have requested for edge " + std::to_string(o.trim) + ".");
            }
            
            x.l.end   -= o.trim;
            x.l.start += o.trim;

            // Name of the region? Locus of the region if not specified.
            const auto rkey = !x.name.empty() ? x.name : x.l.key();
            
            A_ASSERT(!c2d[x.cID].r2d.count(rkey));
            
            // Eg: chr1 0 248956422 chr1
            c2d[x.cID].r2d[rkey] = x;
            
            f(x, p);
        });

        return c2d;
    }
    
    inline BedData readRegions(const Reader &r)
    {
        return readRegions(r, [](const ParserBed::Data &, const ParserProgress &) {});
    }
}

#endif
