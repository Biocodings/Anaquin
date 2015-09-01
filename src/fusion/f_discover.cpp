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
                             "Covered:     %9%";

        options.writer->write((boost::format(summary) % file
                                                      % stats.m.skip
                                                      % stats.m.nq
                                                      % stats.m.nr
                                                      % stats.m.sn()
                                                      % stats.m.sp()
                                                      % "????"
                                                      % "????"
                                                      % stats.covered).str());
        options.writer->close();
    }

    /*
     * Generating sequin statistics
     */

    {
        options.info("Generating sequins statistics");
        options.writer->open("FusionDiscover_quins.stats");
        
        const auto summary = "Summary for dataset: %1% :\n\n"
                             "   Detected: %2% (%3%) sequins\n"
                             "   Undetected: %4% (%5%) sequins\n\n"
                             "#--------------------|   normal   |  fuzzy   |  e_cov |  o_cov\n";

        // Proportion of sequins detected
        const auto detect = std::count_if(stats.h.begin(), stats.h.end(), [&](const std::pair<SequinID, Counts> &p)
        {
            return p.second;
        });

        const auto prop = (detect / static_cast<double>(stats.h.size()));
        
        options.writer->write((boost::format(summary) % file
                                                      % detect
                                                      % prop
                                                      % (stats.h.size() - detect)
                                                      % (1 - prop)
                               ).str());

        const auto format  = "    %1%:       %2%     %3%     %4%    %5%";

        for (const auto &i : stats.h)
        {
            if (i.second)
            {
                options.writer->write((boost::format(format) % i.first
                                                             % "yes"
                                                             % "-"
                                                             % stats.cov.at(i.first).x
                                                             % stats.cov.at(i.first).y
                                       ).str());
            }
            else
            {
                options.writer->write((boost::format(format) % i.first
                                                             % "-"
                                                             % "-"
                                                             % "-"
                                                             % "-"
                                       ).str());
            }
        }

        options.writer->close();
    }

    return stats;
}