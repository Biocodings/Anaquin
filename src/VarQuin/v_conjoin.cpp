#include "tools/tools.hpp"
#include "VarQuin/v_conjoin.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

extern Scripts PlotConjoint();

VConjoint::Stats VConjoint::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;
    const auto regs = r.regs1();

    VConjoint::Stats stats;
    
    for (const auto &i : regs)
    {
        stats.data[i.first];
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

}

static void writeDetected(const FileName &file,
                          const VConjoint::Stats &stats,
                          const VConjoint::Options &o)
{
}

static void writeSummary(const FileName &file,
                         const FileName &seqs,
                         const VConjoint::Stats &stats,
                         const VConjoint::Options &o)
{
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
     * Generating VarConjoint_detected.csv
     */
    
    writeDetected("VarConjoint_detected.csv", stats, o);
    
    /*
     * Generating VarConjoint_linear.R
     */
    
    o.generate("VarConjoint_linear.R");
    o.writer->open("VarConjoint_linear.R");
    //o.writer->write(PlotConjoint("VarConjoint_detected.csv", "data$Depth", "'FP'"));
    o.writer->close();
}
