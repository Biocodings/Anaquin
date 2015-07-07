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
                    for (const auto &i : s.f_r_fusions)
                    {
                        id = i.first;
                        
                        if (i.second.start == start)
                        {
                            return Positive;
                        }
                    }
                }
                else
                {
                    for (const auto &i : s.f_f_fusions)
                    {
                        id = i.first;
                        
                        if (i.second.start == start)
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
    stats.p.m.nr = s.f_f_fusions.size() + s.f_r_fusions.size();

    /*
     * Write out the statistics
     */

    //    AnalyzeReporter::report("rna_abundance.stats", "abundance.R", stats, "FPKM", c, options.writer);
    
    AnalyzeReporter::report("fusion_discover.stats", stats.p, stats.c, options.writer);
    
    return stats;
}