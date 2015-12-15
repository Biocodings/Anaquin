#include "variant/v_discover.hpp"
#include "variant/v_classify.hpp"

using namespace Anaquin;

VDiscover::Stats VDiscover::report(const FileName &file, const Options &o)
{
    VDiscover::Stats stats;
    const auto &r = Standard::instance().r_var;

    o.writer->open("VarDiscover_false.stats");
    
    classify(stats, file, o, [&](const VCFVariant &v, const Variation *match)
    {
        // Empty Implementation
    });

    stats.chrT->m.nr() = r.countVars();

    o.info("Generating summary statistics");

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Experiment:  %2% variants\n"
                         "   Synthetic:   %3% variants\n"
                         "   Reference:   %4% variants\n"
                         "   Detected:    %5% variants\n"
                         "   False-Pos:   %6% variants\n\n"
                         "   Sensitivity: %7%\n"
                         "   Specificity: %8%";

    o.writer->open("VarDiscover_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.chrT->n_expT
                                            % stats.chrT->n_chrT
                                            % r.countVars()
                                            % stats.chrT->m.tp()
                                            % (stats.chrT->n_chrT - stats.chrT->m.tp())
                                            % stats.chrT->m.sn()
                                            % stats.chrT->m.ac()).str());
    o.writer->close();

    /*
     * Generate statistics for each sequin
     */

    o.writer->open("VarDiscover_quins.stats");
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "id" % "detected").str());

    for (const auto &i : stats.chrT->h)
    {
        o.writer->write((boost::format(format) % i.first % stats.chrT->h.at(i.first)).str());
    }
    
    o.writer->close();

    return stats;
}