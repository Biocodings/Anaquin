#include "variant/v_discover.hpp"
#include "variant/v_classify.hpp"

using namespace Anaquin;

VDiscover::Stats VDiscover::analyze(const std::string &file, const Options &o)
{
    VDiscover::Stats stats;
    const auto &r = Standard::instance().r_var;

    classify(stats, file, o, [&](const VCFVariant &v, const Variation *match)
    {
        // Empty Implementation
    });
    
    stats.m.nr = r.countVars();

    o.info("Generating statistics");
    
    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Query: %2% variants\n"
                         "   Filtered: %3% variants (in chrT)\n"
                         "   Detected: %4% variants\n"
                         "   False-Pos: %5% variants\n"
                         "   Reference: %6% variants\n\n"
                         "   Fuzzy: %7%\n"
                         "   Sensitivity:\t%8%\n"
                         "   Specificity:\t%9%"
    ;

    o.writer->open("VarDiscover_false.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % stats.detected
                                            % (stats.n_chrT - stats.detected)
                                            % r.countVars()
                                            % o.fuzzy
                                            % stats.m.sn()
                                            % stats.m.sp()
                     ).str());
    o.writer->close();

    return stats;
}