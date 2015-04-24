#include "expression.hpp"
#include "dna/structural.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Spike;

struct VariantSearch
{
    // Whether a match has been found for the position
    bool matched;
    
    // Whether the zygosity is matched
    bool zy;
    
    // Whether the sequence is matched
    bool seq;
    
    // Whether the alleles are matched
    bool alts;

    const Variation *ref;
};

static bool find(const std::vector<Variation> &vs, const VCFVariant &q, VariantSearch &r)
{
    for (auto &v : vs)
    {
        if (v.pos == q.l.start)
        {
            r.zy   = (v.zy == r.zy);
            r.seq  = (v.r == q.r);
            r.alts = (q.alts.size() == 1 && q.alts.count(v.m));
            r.ref  = &v;
            return (r.matched = true);
        }
    }

    return (r.matched = false);
}

StructuralStats Structural::analyze(const std::string &file, const Options &options)
{
    const auto &r = Standard::instance();

    StructuralStats stats;
    VariantSearch vs;
    
    auto cb = countsForGenes();

    ParserVCF::parse(file, [&](const VCFVariant &v)
    {
        classify(stats, v, [&](const VCFVariant &)
        {
            if (r.d_vars.count(v.id) && vs.zy && vs.seq && vs.alts)
            {
                cb[vs.ref->id]++;
                return true;
            }
            else
            {
                return false;
            }
        });
    });

    stats.mb.fn() = stats.nr - stats.mb.tp();

    /*
     * Calculate the limit of sensitivity
     */
    
    const auto rb = Expression::analyze(cb);
    
    stats.sb = rb.sens(r.r_seqs_gA);

    // Report for the base-level
    AnalyzeReporter::reportClassify("variation_base.stats", stats.dilution(), stats.mb, stats.sb, cb, options.writer);

    return stats;
}