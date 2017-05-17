#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

VCopy::Stats VCopy::analyze(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto &r = Standard::instance().r_var;
 
    VCopy::Stats stats;
    
    std::map<Locus, CopyNumber> l2c;
    std::map<CopyNumber, Locus> c2l;
    
    /*
     * Filter all sequins with CNV equal to the input option
     */
    
    for (auto &i : r.data())
    {
        const auto m = r.match(i.first);
        l2c[m->l] = m->concent();
        c2l[m->concent()] = m->l;
    }

    /*
     * Check calibration statistics for all sequins
     */
    
    // Regions to subsample after trimming
    const auto tRegs = r.regions(true);
    
    // Regions without trimming
    const auto regs = r.regions(false);

    // Check calibration statistics
    auto before = VSample::check(endo, seqs, tRegs, regs, o);

    /*
     * Average coverage on those genomic sequins (ignore others)
     */
    
    Proportion sum = 0;
    unsigned n = 0;
    
    for (const auto &i : before.norms)
    {
        for (const auto &j : i.second)
        {
            if (l2c.at(j.first) == o.copy)
            {
                sum += j.second;
                n++;
            }
        }
    }
    
    /*
     * We have normalization for each genomic sequin, but how much to normalize
     * for others? Let's use the average.
     */
    
    // Average normalization for all genomic sequins
    const auto norm = sum / n;
    
    o.logInfo("Average normalization: " + std::to_string(norm));
    
    /*
     * Adjust non-reference sequins relative to the references
     */
    
    for (auto &i : before.norms)
    {
        for (auto &j : i.second)
        {
            if (l2c.at(j.first) != o.copy)
            {
                j.second = norm;
            }
        }
    }
    
    VSample::Stats s_;
    
    // Perform calibration by subsampling
    const auto after = VSample::sample(seqs, before.norms, s_, regs, tRegs, o);
    
    return stats;
}

void VCopy::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    VCopy::analyze(endo, seqs, o);
}
