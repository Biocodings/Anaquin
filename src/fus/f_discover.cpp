#include "fus/f_discover.hpp"
#include "parsers/parser_fusion.hpp"

using namespace Spike;

FDiscover::Stats FDiscover::analyze(const std::string &file, const Options &options)
{
    FDiscover::Stats stats;
    const auto &s = Standard::instance();
    
    ParserFusion::parse(Reader(file), [&](const ParserFusion::Fusion &f, const ParserProgress &)
    {
        SequinID id;

        if (classify(stats.p.m, f, [&](const ParserFusion::Fusion &)
        {
            if (f.chr_1 == s.id && f.chr_2 == s.id)
            {
                const auto start = f.start_1;
                    
                if (f.dir_1 == Backward && f.dir_2 == Forward)
                {
                    for (const auto &i : s.f_f_exons)
                    {
                        id = i.first;
                        
                        if (i.second.region().start == start)
                        {
                            return Positive;
                        }
                    }
                }
                else
                {
                    for (const auto &i : s.f_f_exons)
                    {
                        id = i.first;

                        if (i.second.region().end == start)
                        {
                            return Positive;
                        }
                    }
                }
            }

            return Negative;
        }))
        {
            // Empty Implementation
        }
    });

    // The references are simply the known fusion points
    stats.p.m.nr = s.f_f_exons.size() + s.f_r_exons.size();

    AnalyzeReporter::report("fusion_discover.stats", stats.p, stats.c, options.writer);
    
    return stats;
}