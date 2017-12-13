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
        std::map<Name, Data> r2d;
    };
    
    struct BedData : public std::map<ChrID, BedChrData>
    {
        inline std::set<Name> names() const
        {
            std::set<Name> x;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.r2d)
                {
                    x.insert(j.first);
                }
            }
            
            return x;
        }

        inline Base len(const ChrID &x, const Name &y = "") const
        {
            Base b = 0;
            
            for (const auto &i : at(x).r2d)
            {
                if (y.empty() || i.second.name == y)
                {
                    b += i.second.l.length();
                }
            }
            
            A_ASSERT(b);
            return b;
        }

        inline Base len() const
        {
            return countMap(*this, [&](const ChrID &x, const BedChrData &)
            {
                return len(x);
            });
        }
        
        inline Counts count() const
        {
            return countMap(*this, [&](const ChrID &, const BedChrData &x)
            {
                return x.r2d.size();
            });
        }
        
        inline DIntervals<> inters(const ChrID &x) const
        {
            DIntervals<> r;
            
            for (const auto &i : at(x).r2d)
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
    
    template <typename F> BedData readRegions(const Reader &r, F f, RegionOptions o = RegionOptions(), const std::set<SequinID> *ex = nullptr)
    {
        BedData c2d;
        
        ParserBed::parse(r, [&](ParserBed::Data &x, const ParserProgress &p)
        {
            if (x.l.length() < 2 * o.trim)
            {
                throw std::runtime_error(x.name + " is too narrow. You have requested for edge " + std::to_string(o.trim) + ".");
            }
            else if (ex && !ex->count(x.name))
            {
                return;
            }
            
            x.l.end   -= o.trim;
            x.l.start += o.trim;

            // Name of the region? Locus of the region if not specified.
            const auto rkey = !x.name.empty() ? x.name : x.l.key();
            
            if (c2d[x.cID].r2d.count(rkey))
            {
                throw std::runtime_error("Duplicate in input reference annoation file. " + x.cID + " " + toString(x.l.start-1) + "-" + toString(x.l.end));
            }
            
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
