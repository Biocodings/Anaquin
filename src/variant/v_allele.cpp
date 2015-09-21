#include "variant/v_allele.hpp"
#include "variant/v_classify.hpp"

using namespace Anaquin;

VAllele::Stats VAllele::report(const std::string &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &r = Standard::instance().r_var;
    
    o.writer->open("VarAllele_false.stats");

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

        // Eg: D_1_12_R_373892_G/A
        const auto id = (boost::format("%1%_%2%_%3%_%4%:") % match->id
                                                           % match->ref
                                                           % match->l.start
                                                           % match->alt).str();

        stats.add(id, known, measured);
    });
 
    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Genome: %2% variants\n"
                         "   Query: %3% variants (in chrT)\n"
                         "   Reference: %4% variants\n"
                         "   Detected: %5% variants\n"
                         "   Correlation:\t%6%\n"
                         "   Slope:\t%7%\n"
                         "   R2:\t%8%\n"
                         "   F-statistic:\t%9%\n"
                         "   P-value:\t%10%\n"
                         "   SSM: %11%, DF: %12%\n"
                         "   SSE: %13%, DF: %14%\n"
                         "   SST: %15%, DF: %16%\n"
    ;
    
    const auto lm = stats.linear(false);
    
    stats.sn = static_cast<double>(stats.detected) / r.countVars();

    o.writer->open("VarAllele_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % r.countVars()
                                            % stats.detected
                                            % lm.r
                                            % lm.m
                                            % lm.r2
                                            % lm.f
                                            % lm.p
                                            % lm.ssm
                                            % lm.ssm_df
                                            % lm.sse
                                            % lm.sse_df
                                            % lm.sst
                                            % lm.sst_df).str());
    o.writer->close();

    AnalyzeReporter::scatter(stats, "", "VarAllele", "Expected allele frequency (proportion)", "Measured allele frequency (proportion)", o.writer);

    return stats;
}