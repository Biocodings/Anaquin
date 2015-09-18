#include "variant/v_allele.hpp"
#include "variant/v_classify.hpp"

using namespace Anaquin;

VAllele::Stats VAllele::report(const std::string &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &r = Standard::instance().r_var;
    
    classify(stats, file, o, [&](const VCFVariant &v, const Variation *match)
    {
        // The known coverage for allele frequnece
        const auto known = r.alleleFreq(Mix_1, match->bID);

        // The measured coverage is the number of base calls aligned and used in variant calling
        const auto measured = static_cast<double>(v.dp_a) / (v.dp_r + v.dp_a);

        /*
         * Plotting the relative allele frequency that is established by differences
         * in the concentration of reference and variant DNA standards.
         */

        stats.add(match->id, known, measured);
    });
 
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
                         "   Sensitivity:\t%8%\n\n"
                         "   Correlation:\t%9%\n"
                         "   Slope:\t%10%\n"
                         "   R2:\t%11%\n"
                         "   Adjusted R2:\t%12%\n"
                         "   F-statistic:\t%13%\n"
                         "   P-value:\t%14%\n"
                         "   SSM: %15%, DF: %16%\n"
                         "   SSE: %17%, DF: %18%\n"
                         "   SST: %19%, DF: %20%\n"
    ;
    
    const auto lm = stats.linear();
    
    stats.sn = static_cast<double>(stats.detected) / r.countVars();

    o.writer->open("VarAllele_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % stats.detected
                                            % (stats.n_chrT - stats.detected)
                                            % r.countVars()
                                            % o.fuzzy
                                            % stats.sn
                                            % lm.r
                                            % lm.m
                                            % lm.r2
                                            % lm.ar2
                                            % lm.f
                                            % lm.p
                                            % lm.ssm
                                            % lm.ssm_df
                                            % lm.sse
                                            % lm.sse_df
                                            % lm.sst
                                            % lm.sst_df).str());
    o.writer->close();

    AnalyzeReporter::scatter(stats, "VarAllele", "Expected allele frequency", "Measured allele frequency", o.writer);

    return stats;
}