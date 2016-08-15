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

        template <typename F> static ReaderBam::Stats stats(const FileName &file,
                                                            const C2Intervals &c2l,
                                                            F f)
        {
            ReaderBam::Stats stats;

            // For each chromosome...
            for (const auto &i : c2l)
            {
                Intervals<> x;
                
                for (const auto &inter : i.second.data())
                {
                    const auto &l = inter.second.l();
                    x.add(Interval(l.key(), l));
                }

                stats.gen[i.first] = x;
                stats.syn[i.first] = x;
                stats.gen[i.first].build();
                stats.syn[i.first].build();
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
                    auto inters  = Standard::isSynthetic(x.cID) ? &stats.syn : &stats.gen;
                    auto matched = (*inters)[x.cID].overlap(x.l);

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