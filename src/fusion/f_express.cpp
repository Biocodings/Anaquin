#include "fusion/f_express.hpp"
#include "parsers/parser_tracking.hpp"

using namespace Anaquin;

FExpress::Stats FExpress::analyze(const std::string &file, const Options &options)
{
    FExpress::Stats stats;
    const auto &s = Standard::instance();

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
    {
        // Don't overflow
        const auto fpkm = std::max(0.05, t.fpkm);
        
        // Try to match by names if possible
        const auto *r = s.f_seqs_A.count(t.geneID) ? &(s.f_seqs_A.at(t.geneID)) : nullptr;

        if (r && t.geneID[0] == 'N')
        {
            if (t.fpkm)
            {
                auto it = std::find(stats.z.begin(), stats.z.end(), t.geneID);

                // Expected FPKM for the sequin
                const auto known = r->abund() * r->length / 1000;
                
                // Measured FPKM for the sequin
                const auto measured = fpkm;
                
                if (it != stats.z.end())
                {
                    auto i = std::distance(stats.z.begin(), it);
                    
                    //stats.x[i] += (r->abund() / r->length);
                    stats.y[i] += (measured);
                }
                else
                {
                    //stats.x.push_back(log2(r->abund()));
                    //stats.y.push_back(log2(fpkm));
                    stats.x.push_back(known);
                    stats.y.push_back(measured);
                    stats.z.push_back(t.geneID);
                }
            }
        }
        else
        {
            std::cout << t.trackID << std::endl;
        }
    });

    for (int i = 0; i < stats.z.size(); i++)
    {
        stats.x[i] = log2(stats.x[i]);
        stats.y[i] = log2(stats.y[i]);
    }

        AnalyzeReporter::linear(stats, "FusionExpression", "FPKM", options.writer);
    
    return stats;
}