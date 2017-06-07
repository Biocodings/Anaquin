#include "VarQuin/v_conjoint.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

extern Scripts PlotConjoint();

VConjoint::Stats VConjoint::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    VConjoint::Stats stats;
    
    for (const auto &i : r.l2Seqs())
    {
        stats.data[i];
    }
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
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
    const auto &r = Standard::instance().r_var;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Sequin" % "Unit" % "Stoichiometry" % "CNV" % "Observed").str());

    for (const auto &i : stats.data)
    {
        const auto seq = noLast(i.first, "_");
        o.writer->write((boost::format(format) % seq
                                               % i.first
                                               % r.concent1(seq)
                                               % r.concent2(i.first)
                                               % i.second.measured).str());
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
    const auto format = "-------VarConjoint Summary Statistics\n\n"
                        "-------VarConjoint Output Results\n\n"
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
                        "       P-value:     %9%\n"
                        "       SSM:         %10%, DF: %11%\n"
                        "       SSE:         %12%, DF: %13%\n"
                        "       SST:         %14%, DF: %15%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % seqs              // 1
                                           % MixRef()          // 2
                                           % r.l1Seqs().size() // 3
                                           % stats.size()      // 4
                                           % ls.m              // 5
                                           % ls.r              // 6
                                           % ls.R2             // 7
                                           % ls.F              // 8
                                           % ls.p              // 9
                                           % ls.SSM            // 10
                                           % ls.SSM_D          // 11
                                           % ls.SSE            // 12
                                           % ls.SSE_D          // 13
                                           % ls.SST            // 14
                                           % ls.SST_D          // 15
            ).str());

    o.writer->close();
}

static void writeConjointR(const FileName &file, const VConjoint::Stats &stats, const VConjoint::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createRLinear("VarConjoint_sequins.csv",
                                           o.work,
                                           "Expected CNV vs Observed Abundance",
                                           "Expected CNV",
                                           "Observed Abundance (log2)",
                                           "log2(data$Stoichiometry * data$CNV)",
                                           "log2(data$Observed)",
                                           "input",
                                           true,
                                           PlotConjoint()));
    o.writer->close();
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

    writeConjointR("VarConjoint_linear.R", stats, o);
}
