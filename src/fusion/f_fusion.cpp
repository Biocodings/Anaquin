#include "fusion/f_fusion.hpp"
#include "parsers/parser_fusion.hpp"

using namespace Anaquin;

FFusion::Stats FFusion::analyze(const std::string &file, const Options &options)
{
    FFusion::Stats stats;
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
            assert(!id.empty());
            
            const auto seq = s.f_seqs_A.at(id);
            
            // Known abundance for the fusion
            const auto known = seq.abund() / seq.length;
            
            // Measured abundance for the fusion
            const auto measured = f.reads;
            
            stats.x.push_back(log2f(known));
            stats.y.push_back(log2f(measured));
            stats.z.push_back(id);
        }
    });

    // The references are simply the known fusion points
    stats.p.m.nr = s.f_f_fusions.size() + s.f_r_fusions.size();

    /*
     * Write out the statistics
     */

    options.writer->open("fusion_coverage.R");
    options.writer->write(RWriter::write(stats.x, stats.y, stats.z, "?", 0.0));
    options.writer->close();
    
    AnalyzeReporter::report("fusion_discover.stats", stats.p, stats.c, options.writer);
    
    return stats;
}