#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

extern FileName BedRef();

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
    
    /*
     * Filter all sequins with CNV equal to the input option
     */
    
    for (auto &i : r.seqs())
    {
        l2c[r.locus(i)] = r.concent1(i);
    }

    /*
     * Check calibration statistics for all sequins
     */
    
    // Regions to subsample after trimming
    const auto tRegs = r.regs1();
    
    // Regions without trimming
    const auto regs = r.regs2();

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
            if (l2c.at(j.first) == o.gen)
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
            if (l2c.at(j.first) != o.gen)
            {
                j.second = norm;
                stats.before.c2v.at(i.first).at(j.first).norm = norm;
            }
        }
    }
    
    // Perform calibration by subsampling
    stats.after = VSample::sample(seqs, stats.before.norms, regs, tRegs, o);

    stats.tBefore = VSample::tBefore(stats.before, stats.after);
    stats.tAfter  = VSample::tAfter (stats.before, stats.after);
    stats.sBefore = VSample::sBefore(stats.before, stats.after);
    stats.sAfter  = VSample::sAfter (stats.before, stats.after);

    stats.afterSeqs = VSample::afterSeqsC(tRegs, stats.before.c2v, o);
    
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
                                                   % r.concent1(j.second.rID)
                                                   % j.second.endo
                                                   % j.second.before
                                                   % j.second.after
                                                   % j.second.norm).str());
        }
    }
    
    o.writer->close();
}

static void generateSummary(const FileName &file,
                            const FileName &endo,
                            const FileName &seqs,
                            const VCopy::Stats &stats,
                            const VCopy::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);

    const auto summary = "-------VarCopy Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Alignment file (genome):  %2%\n"
                         "       Alignment file (sequins): %3%\n\n"
                         "-------Reference regions\n\n"
                         "       Sequin regions: %4% regions\n"
                         "       Method: %5%\n\n"
                         "-------Total alignments (before subsampling)\n\n"
                         "       Genome:    %6%\n"
                         "       Synthetic: %7%\n\n"
                         "-------Total alignments (after subsampling)\n\n"
                         "       Genome:    %8%\n"
                         "       Synthetic: %9%\n\n"
                         "-------Alignments within sampling regions (before subsampling)\n\n"
                         "       Genome:    %10%\n"
                         "       Synthetic: %11%\n\n"
                         "-------Alignments within sampling regions (after subsampling)\n\n"
                         "       Genome:    %12%\n"
                         "       Synthetic: %13%\n\n"
                         "       Normalization: %14% \u00B1 %15%\n\n"
                         "-------Before subsampling (within sampling regions)\n\n"
                         "       Genome coverage (average):    %16%\n"
                         "       Synthetic coverage (average): %17%\n\n"
                         "-------After subsampling (within sampling regions)\n\n"
                         "       Genome coverage (average):    %18%\n"
                         "       Synthetic coverage (average): %19%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()                 // 1
                                            % endo                     // 2
                                            % seqs                     // 3
                                            % r.nRegs()                // 4
                                            % meth2Str(o.meth)         // 5
                                            % stats.tBefore.nEndo      // 6
                                            % stats.tBefore.nSeqs      // 7
                                            % stats.tAfter.nEndo       // 8
                                            % stats.tAfter.nSeqs       // 9
                                            % stats.sBefore.nEndo      // 10
                                            % stats.sBefore.nSeqs      // 11
                                            % stats.sAfter.nEndo       // 12
                                            % stats.sAfter.nSeqs       // 13
                                            % stats.before.normMean()  // 14
                                            % stats.before.normSD()    // 15
                                            % stats.before.meanBEndo() // 16
                                            % stats.before.meanBSeqs() // 17
                                            % stats.before.meanBEndo() // 18
                                            % stats.afterSeqs          // 19
                     ).str());
    o.writer->close();
}

void VCopy::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto stats = VCopy::analyze(endo, seqs, o);

    /*
     * Generating VarCopy_summary.stats
     */
    
    generateSummary("VarCopy_summary.stats", endo, seqs, stats, o);
    
    /*
     * Generating VarCopy_sequins.csv
     */
    
    generateCSV("VarCopy_sequins.csv", stats, o);
}
