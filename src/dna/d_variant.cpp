#include "dna/d_variant.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

VariantStats DVariant::analyze(const std::string &file, const Options &options)
{
    VariantStats stats;
    const auto &s = Standard::instance();

    // Create a counter for each known variant
    auto c = DAnalyzer::counterSequins();

    Variation r;

    ParserVCF::parse(file, [&](const VCFVariant &q, const ParserProgress &)
    {
        if (classify(stats.p.m, q, [&](const VCFVariant &)
        {
            // Can we find the known variant?
            if (!s.d_vars.count(q.l))
            {
                return Negative;
            }

            r = s.d_vars.at(q.l);

           
            
           
            bool f_found = false;
            
            for (const auto &var: s.d_vars)
            {
                if (var.first.contains(r.l))
                {
                    f_found = true;
                    break;
                }
            }
            
            
            
            if (!f_found)
            {
                std::cout << "????" << std::endl;
            }
            
            /*
            if (options.filters.count(tt.id))
            {
                return Ignore;
            }

            c.at(tt.name)++;
*/
            //stats.p_l.m.tp()++;
            
            // Does the alternative allele match with the reference?
            const bool match_a  = (q.a == r.a);
            
            // Does the reference allele match with the reference?
            const bool match_r  = (q.r == r.r);
            
            // Does the genotype match with the reference?
            const bool match_gt = (q.gt == r.gt);
            
            // Does the allele frequency? One would expect to match for a perfect experiment.
            //const bool match_af = (q.af == r.af);
            
            //if (match_gt)           { stats.p_gt.m.tp()++; }
            //if (match_af)           { stats.p_af.m.tp()++; }
            //if (match_a && match_r) { stats.p_al.m.tp()++; }
            
            return (match_a & match_r & match_gt ? Positive : Negative);
        }))
        {
            c[r.id]++;
        }
    });

    //stats.p.m.nr = stats.p_l.m.nr = stats.p_gt.m.nr = stats.p_af.m.nr = stats.p_al.m.nr = s.d_vars.size();
    stats.p.m.nr = s.d_vars.size();

    //assert(stats.p_l.m.tp()  >= stats.p.m.tp());
    //assert(stats.p_gt.m.tp() >= stats.p.m.tp());
    //assert(stats.p_af.m.tp() >= stats.p.m.tp());

    // The structure depends on the mixture
    //const auto seqs = s.r_pair(options.mix);

    /*
     * Calculate for the LOS
     */

    //stats.s = Expression::analyze(c, seqs);
    //AnalyzeReporter::report("dalign_overall.stats", stats, stats.p.m, stats.p.s, c, options.writer);

    return stats;
}