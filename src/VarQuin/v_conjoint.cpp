#include "VarQuin/v_conjoint.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

extern Scripts PlotConjoint();

VConjoint::Stats VConjoint::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto regs = r.regs1();
    
    VConjoint::Stats stats;
    
    for (const auto &i : r.seqs())
    {
        stats.data[i];
    }
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        /*
         * Chromosome is really conjoint sequin name here
         */
        
        if (stats.data.count(x.cID))
        {
            stats.data[x.cID].measured++;
        }
    });

    return stats;
}

static void writeQuins(const FileName &file,
                       const VConjoint::Stats &stats,
                       const VConjoint::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name" % "Sequin" % "Expected" % "Observed" % "Fold").str());

    for (const auto &i : stats)
    {
        o.writer->write((boost::format(format) % i.first
                                               % noLast(i.first, "_")
                                               % i.second.x
                                               % i.second.y
                                               % last(i.first, "_")).str());
    }

    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &seqs,
                         const VConjoint::Stats &stats,
                         const VConjoint::Options &o)
{
    const auto &r = Standard::instance().r_var;
    extern FileName MixRef();

    const auto ls = stats.linear();
    const auto format = "-------VarConjoint Output\n\n"
                        "       Summary for input: %1%\n"
                        "       Mixture file: %2%\n\n"
                        "-------Reference Annotations\n\n"
                        "       Sequins: %3%\n"
                        "-------Detected sequins\n\n"
                        "       Sequins: %4%\n"
                        "-------Linear regression (log2 scale)\n\n"
                        "       Slope:       %5%\n"
                        "       Correlation: %6%\n"
                        "       R2:          %7%\n"
                        "       F-statistic: %8%\n"
                        "       P-value:     %9\n"
                        "       SSM:         %10%, DF: %11%\n"
                        "       SSE:         %12%, DF: %13%\n"
                        "       SST:         %14%, DF: %15%\n";
    
    o.writer->write((boost::format(format) % seqs             // 1
                                           % MixRef()         // 2
                                           % r.regs1().size() // 3
                                           % stats.size()     // 4
                                           % ls.m             // 5
                                           % ls.r             // 6
                                           % ls.R2            // 7
                                           % ls.F             // 8
                                           % ls.p             // 9
                                           % ls.SSM           // 10
                                           % ls.SSM_D         // 11
                                           % ls.SSE           // 12
                                           % ls.SSE_D         // 13
                                           % ls.SST           // 14
                                           % ls.SST_D         // 15
            ).str());
}

void VConjoint::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");

    /*
     * Generating VarConjoint_summary.stats
     */
    
    writeSummary("VarConjoint_summary.stats", file, stats, o);

    /*
     * Generating VarConjoint_sequins.csv
     */
    
    writeQuins("VarConjoint_sequins.csv", stats, o);

    /*
     * Generating VarConjoint_linear.R
     */
    
    o.generate("VarConjoint_linear.R");
    o.writer->open("VarConjoint_linear.R");
    //o.writer->write(PlotConjoint("VarConjoint_detected.csv", "data$Depth", "'FP'"));
    o.writer->close();
}
