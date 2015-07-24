#include "fusion/f_discover.hpp"
#include "parsers/parser_fusion.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::analyze(const std::string &file, const Options &options)
{
    FDiscover::Stats stats;
    const auto &s = Standard::instance();

    options.info("Parsing alignment file");

    ParserFusion::parse(Reader(file), [&](const ParserFusion::Fusion &f, const ParserProgress &p)
    {
        if ((p.i % 1000000) == 0)
        {
            options.wait(std::to_string(p.i));
        }

        options.logInfo((boost::format("%1%: %2% %3%") % p.i % f.chr_1 % f.chr_2).str());
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

            stats.h[seq.id]++;
            stats.x.push_back(log2f(known));
            stats.y.push_back(log2f(measured));
            stats.z.push_back(id);
        }
    });

    // The references are simply the known fusion points
    stats.p.m.nr = s.f_f_fusions.size() + s.f_r_fusions.size();

    options.info("Calculating limit of sensitivity");
    stats.p.s = Expression::analyze(stats.h, s.f_seqs_A);
    stats.s = stats.p.s; // TODO: Fix this

    options.info("Generating statistics");

    AnalyzeReporter::linear(stats, "fusion_align", "FPKM", options.writer);
    AnalyzeReporter::stats("fusion_sequins.stats", stats.p, stats.h, options.writer);

    return stats;
}