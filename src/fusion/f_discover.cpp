#include "fusion/f_discover.hpp"
#include "fusion/f_analyzer.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::analyze(const std::string &file, const FDiscover::Options &options)
{
    const auto stats = FAnalyzer::analyze<FDiscover::Options, FDiscover::Stats>(file, options);

    /*
     * Generate summary statistics
     */

    { 
        options.info("Generating summary statistics");
        options.writer->open("FusionDiscover_summary.stats");

        const auto summary = "Summary for dataset: %1% :\n\n"
                             "   Ignored: %2% fusions not in chrT\n"
                             "   Query: %3% fusions in chrT\n"
                             "Reference: %4% sequins\n\n"
                             "#--------------------|   Sn   |  Sp   |  fSn |  fSp\n"
                             "    Fusion level:       %5%     %6%     %7%    %8%\n"
                             "\n"
                             "Covered:     %9%"
        ;

        options.writer->write((boost::format(summary) % file
                                                      % stats.m.skip
                                                      % stats.m.nq
                                                      % stats.m.nr
                                                      % stats.m.sp()
                                                      % stats.m.sn()
                                                      % "????"
                                                      % "????"
                                                      % stats.covered).str());
        options.writer->close();
    }

    /*
     * Generating sequin statistics
     */

    {
        const auto format = "%1%";

        options.info("Generating sequins statistics");
        options.writer->open("FusionDiscover_quins.stats");
        
        for (const auto &i : stats.h)
        {
            options.writer->write((boost::format(format) % i.first).str());
        }

        options.writer->close();
    }

    {
        AnalyzeReporter::linear(stats, "FusionDiscoverAbund", "FPKM", options.writer);
    }
    
    {
        AnalyzeReporter::missing("FusionDiscover_miss.csv", stats, options.writer);
    }

    return stats;
}