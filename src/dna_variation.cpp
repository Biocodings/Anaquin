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
    
    bool mTrue = false;
    
    unsigned tp = 0;
    unsigned fp = 0;
    unsigned fn = 0;
    unsigned tn = 0;
    
    ParserVCF::parse("tests/data/dna_sims/DNA.flat.chrT.vcf", [&](const VCFVariant &v)
    {
        if (v.chID == r.id)
        {
            if (find(r.vars, v, mTrue))
            {
                n++;
                
                if (mTrue)
                {
                    tp++;
                    nt++;
                }
                else
                {
                    fp++;
                }
            }
            else
            {
                n++;
                fp++;
            }
        }
        else
        {
            if (find(r.vars, v, mTrue))
            {
                //if (mTrue)
                {
                    fn++;
                }
            }
            else
            {
                tn++;
            }
        }
    });

    VariationStats stats;

    stats.covered = static_cast<Percentage>(n) / r.vars.size();
    stats.efficiency = static_cast<Percentage>(nt) / n;
    
    return stats;
}