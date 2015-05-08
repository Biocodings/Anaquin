#include "expression.hpp"
#include "dna/structural.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

StructuralStats Structural::analyze(const std::string &file, const Options &options)
{
    StructuralStats stats;
    const auto &s = Standard::instance();
    
    auto c = countsForGenes();
    Variation r;

    ParserVCF::parse(file, [&](const VCFVariant &q, const ParserProgress &)
    {
        if (q.m == SNP || q.m == Indel)
        {
            /*
             * SNP/Indel classification
             */

            if (classify(stats.m, q, [&](const VCFVariant &)
            {
                if (!s.d_vars.count(q.l))
                {
                    return Negative;
                }

                // Can we find a matching reference variation?
                r = s.d_vars.at(q.l);

                if (options.filters.count(r.id))
                {
                    return Ignore;
                }

                return (q.r == r.r && q.a == r.a ? Positive : Negative);
            }))
            {
                c[r.id]++;
            }
        }
    });

    
    
    
    
    
    
    
    // The structure depends on the mixture
    const auto seqs = s.r_pair(options.mix);

    
    
    /*
     * Calculate for the LOS
     */
    
    stats.s = Expression::analyze(c, seqs);

    
    
    
    //stats.mb.fn() = stats.nr - stats.mb.tp();

    /*
     * Calculate the limit of sensitivity
     */
    
//    const auto rb = Expression::analyze(cb);
//    
//    stats.sb = rb.sens(r.r_seqs_gA);
//
//    // Report for the base-level
//    AnalyzeReporter::report("variation_base.stats", stats.dilution(), stats.mb, stats.sb, cb, options.writer);
//
    return stats;
}