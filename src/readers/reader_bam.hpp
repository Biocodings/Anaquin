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
            ID2Intervals inters;
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

                stats.inters[i.first] = x;
                stats.inters[i.first].build();
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
                
                if (x.mapped && stats.inters.count(x.cID))
                {
                    auto matched = stats.inters[x.cID].overlap(x.l);

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