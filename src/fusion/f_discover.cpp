#include "fusion/f_discover.hpp"
#include "fusion/f_classify.hpp"

using namespace Anaquin;

FDiscover::Stats FDiscover::report(const std::string &file, const FDiscover::Options &o)
{
    const auto stats = FClassify::analyze<FDiscover::Options, FDiscover::Stats>(file, false, o);

    /*
     * Generate summary statistics
     */

    { 
        o.info("Generating summary statistics");
        o.writer->open("FusionDiscover_summary.stats");

        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Genome: %2% fusions\n"
                             "   Synthetic: %3% fusions\n"
                             "   Genome-Synthetic: %4% fusions\n"
                             "   Reference: %5% sequins\n\n"
                             "   Fuzzy: %6%\n\n"
                             "   Sensitivity: %7%\n"
                             "   Specificity: %8%\n";

        o.writer->write((boost::format(summary) % file
                                                % stats.n_hg38
                                                % stats.n_chrT
                                                % stats.hg38_chrT
                                                % stats.m.nr
                                                % o.fuzzy
                                                % stats.m.sn()
                                                % stats.m.sp()).str());
        o.writer->close();
    }

    /*
     * Generating statistics for sequins
     
     */

    {
        o.info("Generating sequins statistics");
        o.writer->open("FusionDiscover_quins.stats");

        const auto format = "%1%\t%2%";

        o.writer->write((boost::format(format) % "id" % "counts").str());

        for (const auto &i : stats.h)
        {
                o.writer->write((boost::format(format) % i.first
                                                       % i.second
                                       ).str());
        }

        o.writer->close();
    }

    return stats;
}