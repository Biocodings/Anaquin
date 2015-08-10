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
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

        options.info("Generate summary statistics");
        options.writer->open("FusionDiscover_summary.stats");
        options.writer->write((boost::format(format) % "sn"
                                                     % "sp"
                                                     % "coverage"
                                                     % "los"
                                                     % "partner_5"
                                                     % "partner_3").str());
        options.writer->write((boost::format(format) % stats.m.sn()
                                                     % stats.m.sp()
                                                     % stats.covered
                                                     % stats.s.abund
                                                     % "?"
                                                     % "?").str());
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
        AnalyzeReporter::linear(stats, "FusionDiscover", "FPKM", options.writer);
    }
    
    {
        AnalyzeReporter::missing("FusionDiscover_miss.csv", stats, options.writer);
    }

    return stats;
}