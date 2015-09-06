#include "variant/v_allele.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

VAllele::Stats VAllele::analyze(const std::string &file, const Options &o)
{
    VAllele::Stats stats;
    const auto &r = Standard::instance().r_var;

    long queries  = 0;
    long filtered = 0;
    long detected = 0;

    o.info("Parsing VCF file");
    o.writer->open("VarAllele_false.stats");
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.writer->write((boost::format(format) % "start"
                                           % "matched"
                                           % "type"
                                           % "alt"
                                           % "ref").str());

    ParserVCF::parse(file, [&](const VCFVariant &v, const ParserProgress &)
    {
        queries++;

        if (v.id != Standard::instance().id)
        {
            return;
        }
        
        filtered++;
        const Variation *match;

        Confusion m;

        if (classify(m, v, [&](const VCFVariant &)
        {
            const auto found = (match = r.findVar(v.l)) != nullptr;
            const auto type  = (match && match->type == v.type);
            const auto alt   = (match && match->alt  == v.alt);
            const auto ref   = (match && match->ref  == v.ref);

            if (!found || !type || !alt || !ref)
            {
                o.writer->write((boost::format(format) % v.l.start
                                                       % found
                                                       % type
                                                       % alt
                                                       % ref).str());
                return Negative;
            }
            
            detected++;
            
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
            stats.h.at(match->id)++;
        }
    });

    o.writer->close();
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
                         "Fuzzy: %7%\n"
                         "Sensitivity:\t%8%\n\n"
                         "Correlation:\t%9%\n"
                         "Slope:\t%10%\n"
                         "R2:\t%11%\n"
                         "Adjusted R2:\t%12%\n"
                         "F-statistic:\t%13%\n"
                         "P-value:\t%14%\n"
                         "SSM: %15%, DF: %16%\n"
                         "SSE: %17%, DF: %18%\n"
                         "SST: %19%, DF: %20%\n"
    ;
    
    const auto lm = stats.linear();
    
    stats.sn = static_cast<double>(detected) / r.countVars();

    o.writer->open("VarAllele_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % queries
                                            % filtered
                                            % detected
                                            % (filtered - detected)
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

    AnalyzeReporter::scatter(stats, "VarAllele", "", o.writer);

    return stats;
}