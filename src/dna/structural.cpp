#include "expression.hpp"
#include "r_analyzer.hpp"
#include "dna/structural.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

StructuralStats Structural::analyze(const std::string &file, const Options &options)
{
    StructuralStats stats;
    const auto &s = Standard::instance();

    auto c = DAnalyzer::d_sequinCounter();
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
                
                stats.ml.tp()++;

                // Can we find a matching reference variation?
                r = s.d_vars.at(q.l);

                if (options.filters.count(r.id))
                {
                    return Ignore;
                }

                return ((q.r == r.r && q.a == r.a) ? Positive : Negative);
            }))
            {
                c[r.id]++;
            }
        }
    });

    stats.m.nr = stats.ml.nr = s.d_vars.size();

    assert(stats.ml.tp() >= stats.m.tp());

    // The structure depends on the mixture
    //const auto seqs = s.r_pair(options.mix);

    /*
     * Calculate for the LOS
     */

    //stats.s = Expression::analyze(c, seqs);

    AnalyzeReporter::report("dalign_overall.stats", stats.m, stats.s, c, options.writer);

    return stats;
}