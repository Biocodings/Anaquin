#include "parser_vcf.hpp"
#include "dna_variation.hpp"
#include "standard_factory.hpp"

static bool find(const std::vector<Variation> &vs, VCFVariant x, bool &mTrue)
{
    mTrue = false;
    
    for (auto v : vs)
    {
        if (v.l.end == x.pos)
        {
            const auto x_r = x.ref;
            const auto x_q = x.alts[0]; // Assume only one for now...
            
            mTrue = (v.r == x_r && v.m == x_q);
            
            return true;
        }
    }
    
    return false;
}

VariationStats DNAVariation::analyze(const std::string &file)
{
    const auto r = StandardFactory::reference();

    unsigned n = 0;
    unsigned nt = 0;
    
    bool mTrue;
    
    ParserVCF::parse("tests/data/simulation/DNA.flat.chrT.vcf", [&](const VCFVariant &v)
    {
        if (find(r.vars, v, mTrue))
        {
            n++;
            
            if (mTrue)
            {
                nt++;
            }
        }
        else
        {
            
        }
    });

    VariationStats stats;

    stats.covered = static_cast<Percentage>(n) / r.vars.size();
    stats.efficiency = static_cast<Percentage>(nt) / n;
    
    return stats;
}