#include "variant/v_allele.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

VAllele::Stats VAllele::analyze(const std::string &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &r = Standard::instance().r_var;

    long n = 0;

    o.info("Parsing VCF file");

    ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
    {
        n++;

        if (v.id != Standard::instance().id)
        {
            return;
        }
        
        const Variation *match;

        Confusion m;

        if (classify(m, v, [&](const VCFVariant &)
        {
            // Can we find this variant?
            if (!(match = r.findVar(v.l)))
            {
                return Negative;
            }

            // Does the variant match with the meta?
            else if (match->type != v.type || match->alt != v.alt || match->ref != v.ref)
            {
                return Negative;
            }

            // The known coverage for allele frequnece
            const auto known = r.alleleFreq(MixA, match->bID);

            // The measured coverage is the number of base calls aligned and used in variant calling
            const auto measured = (double) v.dp_a / (v.dp_r + v.dp_a);
            
            /*
             * Plotting the relative allele frequency that is established by differences
             * in the concentration of reference and variant DNA standards.
             */

            stats.add(match->id, known, measured);
  
            return Positive;
        }))
        {
            //stats.h.at(match)++;
        }
    });

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    const auto summary = "Summary for dataset: %1% :\n\n"
                         "   Ignored: %2% variants not in chrT\n"
                         "   Detected: %3% variants in chrT\n"
                         "   Reference: %4% variants\n\n"
                         "Fuzzy: %5%\n\n"
                         "Correlation:     %6%\n"
                         "Slope:     %7%\n"
                         "R2:     %8%\n"
    ;
    
    const auto lm = stats.linear();
    
    o.writer->open("VarAllele_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % (n - stats.size())
                                            % stats.size()
                                            % r.countVars()
                                            % o.fuzzy
                                            % lm.r
                                            % lm.m
                                            % lm.r2).str());
    o.writer->close();

    AnalyzeReporter::scatter(stats, "VarAllele", "", o.writer);

    return stats;
}