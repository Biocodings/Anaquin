#include "fusion/f_discover.hpp"
#include "fusion/f_analyzer.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::analyze(const std::string &file, const FDiscover::Options &options)
{
    const auto stats = FAnalyzer::analyze(file, options);

    /*
     * Generating summary statistics
     */

    {
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";
        
        options.info("Generating summary statistics");
        options.writer->open("FusionDiscover_summary.stats");
        options.writer->write((boost::format(format) % "sn"
                               % "sp"
                               % "coverage"
                               % "los"
                               % "partner_5"
                               % "partner_3").str());
        options.writer->write((boost::format(format) % stats.m.sn()
                               % stats.m.sp()
                               % stats.m.sp()
                               % stats.s.abund
                               % 0
                               % 0).str());
        options.writer->close();
    }
    
    /*
     * Generating sequin statistics
     */

    {
        const auto format = "%1%";

        options.info("Generating sequins statistics");
        options.writer->open("FusionDiscover_quins.stats");
        
        std::cout << stats.h.size() << std::endl;
        
        for (const auto &i : stats.h)
        {
            options.writer->write((boost::format(format) % "id").str());
            options.writer->write((boost::format(format) % i.first).str());
        }
       
        options.writer->close();
    }

    /*
     * The statistics we need are the subset that we've got, excluding the linear modeling
     */
    
    FDiscover::Stats d_stats;
    
    d_stats.m = stats.m;
    d_stats.h = stats.h;
    
    return d_stats;
}