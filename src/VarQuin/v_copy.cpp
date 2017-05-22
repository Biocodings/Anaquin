#include "tools/tools.hpp"
#include "VarQuin/v_copy.hpp"

using namespace Anaquin;

extern FileName BedRef();
extern Scripts PlotCNV();

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

    const auto r1 = r.regs1();
    const auto r2 = r.regs2();

    // Check calibration statistics
    stats.before = VSample::check(endo, seqs, r2, r1, o);

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
     * We have normalization for the sequins, but how much to normalize for others?
     */
    
    // Average normalization for other sequins
    stats.gNorm = sum / n;
    
    o.logInfo("Average normalization: " + std::to_string(stats.gNorm));
    
    /*
     * Adjust non-reference sequins relative to the references
     */
    
    for (auto &i : stats.before.norms)
    {
        for (auto &j : i.second)
        {
            if (l2c.at(j.first) != o.gen)
            {
                j.second = stats.gNorm;
                stats.before.c2v.at(i.first).at(j.first).norm = stats.gNorm;
            }
        }
    }
    
    // Perform calibration by subsampling
    stats.after = VSample::sample(seqs, stats.before.norms, r1, r2, o);
    
    stats.tBefore = VSample::tBefore(stats.before, stats.after);
    stats.tAfter  = VSample::tAfter (stats.before, stats.after);
    stats.sBefore = VSample::sBefore(stats.before, stats.after);
    stats.sAfter  = VSample::sAfter (stats.before, stats.after);

    stats.afterSeqs = VSample::afterSeqsC(r2, stats.before.c2v, o);

    /*
     * Quantifying the CNV ladder
     */
    
    for (const auto &i : countForID(r2))
    {
        const auto exp = r.concent1(i.first);
        const auto obs = i.second;
        stats.add(i.first, exp, log2(obs));
    }
    
    return stats;
}

static void writeCSV(const FileName &file, const VCopy::Stats &stats, const VSample::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);
    
    const auto format = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%");
    
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
                                           % "CNV"
                                           % "Observed"
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
                                                   % stats.at(j.second.rID).x
                                                   % stats.at(j.second.rID).y
                                                   % j.second.norm).str());
        }
    }
    
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &endo,
                         const FileName &seqs,
                         const VCopy::Stats &stats,
                         const VCopy::Options &o)
{
    const auto &r = Standard::instance().r_var;

    o.generate(file);

    const auto lm = stats.linear(false);
    const auto summary = "-------VarCopy Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Alignment file (sample):   %2%\n"
                         "       Alignment file (sequin):   %3%\n\n"
                         "-------Reference regions\n\n"
                         "       Sequin regions: %4% regions\n"
                         "       Method: %5%\n\n"
                         "-------CNV Normalization\n\n"
                         "       Genomic CNV: 2x\n"
                         "       Genomic normalization: %6%\n\n"
                         "-------Overall linear regression\n\n"
                         "       Slope:       %7%\n"
                         "       Correlation: %8%\n"
                         "       R2:          %9%\n"
                         "       F-statistic: %10%\n"
                         "       P-value:     %11%\n"
                         "       SSM:         %12%, DF: %13%\n"
                         "       SSE:         %14%, DF: %15%\n"
                         "       SST:         %16%, DF: %17%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()         // 1
                                            % endo             // 2
                                            % seqs             // 3
                                            % r.nRegs()        // 4
                                            % meth2Str(o.meth) // 5
                                            % stats.gNorm      // 6
                                            % lm.m             // 7
                                            % lm.r             // 8
                                            % lm.R2            // 9
                                            % lm.F             // 10
                                            % lm.p             // 11
                                            % lm.SSM           // 12
                                            % lm.SSM_D         // 13
                                            % lm.SSE           // 14
                                            % lm.SSE_D         // 15
                                            % lm.SST           // 16
                                            % lm.SST_D         // 17
                     ).str());
    o.writer->close();
}

void VCopy::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto stats = VCopy::analyze(endo, seqs, o);

    /*
     * Generating VarCopy_summary.stats
     */
    
    writeSummary("VarCopy_summary.stats", endo, seqs, stats, o);
    
    /*
     * Generating VarCopy_sequins.csv
     */
    
    writeCSV("VarCopy_sequins.csv", stats, o);

    /*
     * Generating VarCopy_linear.R
     */

    o.generate("VarCopy_linear.R");
    o.writer->open("VarCopy_linear.R");
    o.writer->write(RWriter::createScript("VarCopy_sequins.csv", PlotCNV()));
    o.writer->close();
}
