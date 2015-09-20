#include "fusion/f_discover.hpp"
#include "fusion/f_classify.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::report(const std::string &file, const FDiscover::Options &o)
{
    const auto stats = FClassify::analyze<FDiscover::Options, FDiscover::Stats>(file, o);

    /*
     * Generate summary statistics
     */

    { 
        o.info("Generating summary statistics");
        o.writer->open("FusionDiscover_summary.stats");

        const auto summary = "Summary for dataset: %1% :\n\n"
                             "   Genome: %2% fusions\n"
                             "   Synthetic: %3% fusions\n"
                             "   Reference: %4% sequins\n\n"
                             "   Fuzzy: %5%\n\n"
                             "#--------------------|   Sn   |  Sp  \n"
                             "    Fusion level:       %6%     %7%  \n"
                             "\n"
                             "Sensitivity:     %8%";

        o.writer->write((boost::format(summary) % file
                                                % stats.m.skip
                                                % stats.m.nq
                                                % stats.m.nr
                                                % o.fuzzy
                                                % stats.m.sn()
                                                % stats.m.sp()
                                                % stats.covered).str());
        o.writer->close();
    }

    /*
     * Generating sequin statistics
     */

    {

        o.info("Generating sequins statistics");
        o.writer->open("FusionDiscover_quins.stats");
        
        const auto summary = "Summary for dataset: %1% :\n\n"
                             "   Detected: %2% (%3%) sequins\n"
                             "   Undetected: %4% (%5%) sequins\n\n"
                             "   Fuzzy: %6%\n\n"
                             "#--------------------|   normal   |  fuzzy   |  expect (attomol/ul) |  measure (reads)\n";

        // Proportion of sequins detected
        const auto detect = std::count_if(stats.h.begin(), stats.h.end(), [&](const std::pair<SequinID, Counts> &p)
        {
            return p.second;
        });

        const auto prop = (detect / static_cast<double>(stats.h.size()));
        
        o.writer->write((boost::format(summary) % file
                                                % detect
                                                % prop
                                                % (stats.h.size() - detect)
                                                % (1 - prop)
                                                % o.fuzzy
                               ).str());

        const auto format  = "    %1%:       %2%     %3%     %4%    %5%";

        for (const auto &i : stats.h)
        {
            if (i.second)
            {
                o.writer->write((boost::format(format) % i.first
                                                       % "yes"
                                                       % "-"
                                                       % stats.at(i.first).x
                                                       % stats.at(i.first).y
                                       ).str());
            }
            else
            {
                o.writer->write((boost::format(format) % i.first
                                                       % "-"
                                                       % "-"
                                                       % "-"
                                                       % "-"
                                       ).str());
            }
        }

        o.writer->close();
    }

    return stats;
}