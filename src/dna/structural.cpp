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
        if (v.pos == q.pos)
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

StructuralStats Structural::analyze(const std::string &file, const Structural::Options &options)
{
    const auto &r = Standard::instance();

    StructuralStats stats;
    VariantSearch vs;

    ParserVCF::parse(file, [&](const VCFVariant &v)
    {
        stats.n++;

        if (v.id == r.id)
        {
            stats.nr++;

            if (find(r.vars, v, vs) && vs.zy && vs.seq && vs.alts)
            {
                stats.m.tp++;
            }
            else
            {
                stats.m.fp++;
            }
        }
        else
        {
            stats.nq++;
        }
    });

    stats.m.fn = stats.nr - stats.m.tp;
    
    return stats;
}