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
        l2c[r.locus(i)] = r.input1(i);
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
    
    /*
     * Genomic coverage for sequins and genome
     */

    std::set<double> cEndo, bSeqs, aSeqs;
    VSample::afterSeqsC(r2, stats.before.c2v, o);
    
    for (const auto &i : stats.before.c2v)
    {
        for (const auto &j : i.second)
        {
            // Genomic region?
            if (l2c.at(j.first) == o.gen)
            {
                cEndo.insert(j.second.endo);
                bSeqs.insert(j.second.before);
                aSeqs.insert(j.second.after);
            }
        }
    }
    
    stats.gEndo = mean(cEndo);
    stats.bSeqs = mean(bSeqs);
    stats.aSeqs = mean(aSeqs);

    /*
     * Quantifying the CNV ladder
     */

    for (const auto &i : stats.before.c2v)
    {
        for (const auto &j : i.second)
        {
            const auto exp = r.input1(j.second.rID);
            const auto obs = i.second;
            stats.add(j.second.rID, exp, j.second.after);
        }
    }

    return stats;
}

static void writeQuins(const FileName &file, const VCopy::Stats &stats, const VSample::Options &o)
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
                                                   % r.input1(j.second.rID)
                                                   % j.second.endo
                                                   % j.second.before
                                                   % j.second.after
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
                         "       Reference annotation file:              %1%\n"
                         "       Alignments (sample):                    %2%\n"
                         "       Alignments (sequin):                    %3%\n\n"
                         "-------Reference regions\n\n"
                         "       Reference regions:                      %4% regions\n"
                         "       Coverage of sequins (2n):               %5%\n"
                         "       Coverage of counterpart genome regions: %6%\n"
                         "       Method:                                 %7%\n"
                         "       Scaling factor:                         %8%\n\n"
                         "-------Overall linear regression\n\n"
                         "       Slope:       %9%\n"
                         "       Correlation: %10%\n"
                         "       R2:          %11%\n"
                         "       F-statistic: %12%\n"
                         "       P-value:     %13%\n"
                         "       SSM:         %14%, DF: %15%\n"
                         "       SSE:         %16%, DF: %17%\n"
                         "       SST:         %18%, DF: %19%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()         // 1
                                            % endo             // 2
                                            % seqs             // 3
                                            % r.nRegs()        // 4
                                            % stats.aSeqs      // 5
                                            % stats.gEndo      // 6
                                            % stats.gNorm      // 7
                                            % meth2Str(o.meth) // 8
                                            % lm.m             // 9
                                            % lm.r             // 10
                                            % lm.R2            // 11
                                            % lm.F             // 12
                                            % lm.p             // 13
                                            % lm.SSM           // 14
                                            % lm.SSM_D         // 15
                                            % lm.SSE           // 16
                                            % lm.SSE_D         // 17
                                            % lm.SST           // 18
                                            % lm.SST_D         // 19
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
    
    writeQuins("VarCopy_sequins.csv", stats, o);

    /*
     * Generating VarCopy_linear.R
     */

    o.generate("VarCopy_linear.R");
    o.writer->open("VarCopy_linear.R");
    o.writer->write(RWriter::createScript("VarCopy_sequins.csv", PlotCNV()));
    o.writer->close();
}
