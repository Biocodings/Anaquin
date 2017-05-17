#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

template <typename Stats> Coverage stats2cov(const VSample::Method meth, const Stats &stats)
{
    switch (meth)
    {
        case VSample::Method::Mean:   { return stats.mean; }
        case VSample::Method::Median: { return stats.p50;  }
        case VSample::Method::Prop:
        case VSample::Method::Reads: { return stats.mean; }
    }
}

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
    stats.before = VSample::check(endo, seqs, tRegs, regs, o);

    /*
     * Average coverage on those genomic sequins (ignore others)
     */
    
    Proportion sum = 0;
    unsigned n = 0;
    
    for (const auto &i : stats.before.norms)
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
    
    for (auto &i : stats.before.norms)
    {
        for (auto &j : i.second)
        {
            if (l2c.at(j.first) != o.copy)
            {
                j.second = norm;
                stats.before.c2v.at(i.first).at(j.first).norm = norm;
            }
        }
    }
    
    VSample::Stats s_;
    
    // Perform calibration by subsampling
    stats.after = VSample::sample(seqs, stats.before.norms, s_, regs, tRegs, o);
    
    std::vector<double> allAfterSeqsC;
    
    /*
     * Assume our subsampling is working, let's check coverage for the regions.
     */
    
    // For each chromosome...
    for (auto &i : tRegs)
    {
        // For each region...
        for (auto &j : i.second.data())
        {
            //stats.count++;
            
            // Coverage after subsampling
            const auto cov = stats2cov(o.meth, j.second.stats());
            
            stats.before.c2v[i.first].at(j.second.l()).after = cov; // ????
            
            // Required for generating summary statistics
            allAfterSeqsC.push_back(cov);
        }
    }

    return stats;
}

static void generateCSV(const FileName &file, const VCopy::Stats &stats, const VSample::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);
    
    const auto format = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%");
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Start"
                                           % "End"
                                           % "Copy"
                                           % "Genome"
                                           % "Before"
                                           % "After"
                                           % "Norm").str());
    
    // For each chromosome...
    for (const auto &i : stats.before.c2v)
    {
        // For each variant...
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % j.second.rID
                                                   % i.first
                                                   % j.first.start
                                                   % j.first.end
                                                   % r.match(j.second.rID)->concent()
                                                   % j.second.endo
                                                   % j.second.before
                                                   % j.second.after
                                                   % j.second.norm).str());
        }
    }
    
    o.writer->close();
}

void VCopy::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto stats = VCopy::analyze(endo, seqs, o);
    
    /*
     * Generating VarCopy_sequins.csv
     */
    
    generateCSV("VarCopy_sequins.csv", stats, o);
}
