#include "tools/tools.hpp"
#include "VarQuin/v_conjoint.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

extern Scripts PlotConjoint();

VConjoint::Stats VConjoint::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    
    VConjoint::Stats stats;
    
    for (const auto &i : r.seqsL2())
    {
        stats.data[i];
    }
    
    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &)
    {
        if (stats.data.count(x.cID))
        {
            stats.data[x.cID]++;
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
    o.writer->write((boost::format(format) % "Name" % "Unit" % "Stoichiometry" % "CNV" % "Observed").str());

    for (const auto &i : stats.data)
    {
        o.writer->write((boost::format(format) % r.t2(i.first)
                                               % i.first
                                               % r.input1(r.t2(i.first))
                                               % r.input2(i.first)
                                               % i.second).str());
    }

    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &seqs,
                         const VConjoint::Stats &stats,
                         const VConjoint::Options &o)
{
    const auto &r = Standard::instance().r_var;
    extern FileName ConRef();

    const auto format = "-------VarConjoint Summary Statistics\n\n"
                        "-------VarConjoint Output Results\n\n"
                        "       Summary for input: %1%\n"
                        "       Mixture file: %2%\n\n"
                        "-------Reference Annotations\n\n"
                        "       Sequins: %3%\n\n"
                        "-------Detected sequins\n\n"
                        "       Sequins: %4%\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % seqs                // 1
                                           % ConRef()            // 2
                                           % r.seqsL1().size()   // 3
                                           % nonZero(stats.data) // 4
            ).str());
    o.writer->close();
}

static void writeConjointR(const FileName &file, const VConjoint::Stats &stats, const VConjoint::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createRConjoint("VarConjoint_sequins.csv",
                                             PlotConjoint(),
                                             o.work,
                                             "Expected CNV vs Observed Abundance",
                                             "Expected CNV",
                                             "Observed Abundance (log2)",
                                             "log2(data$Stoichiometry * data$CNV)",
                                             "log2(data$Observed)"));
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
