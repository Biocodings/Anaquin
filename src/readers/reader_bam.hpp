#ifndef READER_BAM_HPP
#define READER_BAM_HPP

#include "data/intervals.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"

namespace Anaquin
{
    struct ReaderBam
    {
        struct Stats : public AlignmentStats
        {
            ID2Intervals inters;
        };

        enum class Response
        {
            OK,
            SKIP_MATCH,
            SKIP_EVERYTHING
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
                Interval *matched = nullptr;
                
                if (x.mapped && stats.inters.count(x.cID))
                {
                    matched = stats.inters[x.cID].overlap(x.l);
                }
                
                const auto r = f(x, info, matched);
                
                if (r != Response::SKIP_EVERYTHING)
                {
                    if (r != Response::SKIP_MATCH && matched)
                    {
                        matched->map(x.l);
                    }
                    
                    if (isVarQuin(x.cID))
                    {
                        stats.nSyn++;
                    }
                    else if (x.cID != "*")
                    {
                        stats.nGen++;
                    }
                    else
                    {
                        stats.nNA++;
                    }
                }
            });
            
            return stats;
        }
    };
}

#endif
