#ifndef READER_BAM_HPP
#define READER_BAM_HPP

#include "data/minters.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct ReaderBam
    {
        struct Stats : public AlignmentStats
        {
            MC2Intervals inters;
        };

        template <typename F> static ReaderBam::Stats stats(const FileName &file,
                                                            const std::map<ChrID, std::set<Locus>> &c2l,
                                                            F f)
        {
            ReaderBam::Stats stats;

            for (const auto &inters : c2l)
            {
                MergedIntervals<> mi;
                
                for (const auto &inter : inters.second)
                {
                    mi.add(MergedInterval(std::to_string(inter.start) + "_" + std::to_string(inter.end), inter));
                }

                stats.inters[inters.first] = mi;
            }

            ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
            {
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
                
                // Eg: track progress
                f(x, info);
            });
            
            return stats;
        }
    };
}

#endif