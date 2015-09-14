#include "fusion/f_express.hpp"
#include "parsers/parser_tracking.hpp"

using namespace Anaquin;

FExpress::Stats FExpress::report(const std::string &file, const Options &options)
{
    FExpress::Stats stats;
//    const auto &s = Standard::instance();
//
//    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &p)
//    {
//        // Don't overflow
//        const auto fpkm = std::max(0.05, t.fpkm);
//
//        // Try to match by names if possible
//        const auto r = nullptr; // s.seqs_1.count(t.geneID) ? &(s.seqs_1.at(t.geneID)) : nullptr;
//
//        if (r)
//        {
//            if (t.fpkm)
//            {
//                // Expected FPKM for the sequin
//                const auto known = r->abund() / r->length;
//                
//                // Measured FPKM for the sequin
//                const auto measured = fpkm;
//                
//                if (stats.count(t.geneID))
//                {
//                    stats.at(t.geneID).y += measured;
//                }
//                else
//                {
//                    stats.add(t.geneID, known, measured);
//                }
//            }
//        }
//        else
//        {
//            std::cout << t.trackID << std::endl;
//        }
//    });
//
//    for (auto &i : stats)
//    {
//        i.second.x = log2(i.second.x);
//        i.second.y = log2(i.second.y);
//    }
//
//    AnalyzeReporter::linear(stats, "FusionExpression", "FPKM", options.writer);
//    
    return stats;
}