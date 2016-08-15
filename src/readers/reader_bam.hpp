#ifndef READER_BAM_HPP
#define READER_BAM_HPP

#include "data/intervals.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct ReaderBam
    {
        struct Stats : public AlignmentStats
        {
            // Intervals for genomic reads
            ID2Intervals gen;

            // Intervals for synthetic reads
            ID2Intervals syn;
        };

        static ID2Intervals reg2Inters(const std::map<ChrID, std::set<Locus>> &c2l)
        {
            ID2Intervals r;

            for (const auto &chr : c2l)
            {
                Intervals<> x;
                
                for (const auto &inter : chr.second)
                {
                    x.add(Interval(std::to_string(inter.start) + "_" + std::to_string(inter.end), inter));
                }
                
                r.add(chr.first, x);
            }
            
            return r;
        }

        template <typename F> static ReaderBam::Stats stats(const FileName &file,
                                                            const std::map<ChrID, std::set<Locus>> &c2l,
                                                            F f)
        {
            ReaderBam::Stats stats;

            for (const auto &chr : c2l)
            {
                Intervals<> x;
                
                for (const auto &inter : chr.second)
                {
                    x.add(Interval(std::to_string(inter.start) + "_" + std::to_string(inter.end), inter));
                }

                stats.gen[chr.first] = x;
                stats.syn[chr.first] = x;
            }

            ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
            {
                if (!f(x, info))
                {
                    return;
                }
                
                if (Standard::isSynthetic(x.cID))
                {
                    stats.countSyn++;
                }
                else if (x.cID != "*")
                {
                    stats.countGen++;
                }
                else
                {
                    stats.countNA++;
                }
                
                if (x.mapped && stats.gen.count(x.cID))
                {
                    auto inters  = Standard::isSynthetic(x.cID) ? stats.syn : stats.gen;
                    auto matched = inters[x.cID].overlap(x.l);
                    
                    if (matched)
                    {
                        matched->map(x.l);
                    }
                }
            });
            
            return stats;
        }
    };
}

#endif