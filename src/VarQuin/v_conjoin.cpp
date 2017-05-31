#include "tools/tools.hpp"
#include "VarQuin/v_conjoin.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

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
     * Generating VarConjoin_summary.stats
     */
    
    writeSummary("VarConjoin_summary.stats", file, stats, o);

    /*
     * Generating VarConjoin_sequins.csv
     */
    
    writeQuins("VarConjoin_sequins.csv", stats, o);

    /*
     * Generating VarConjoin_detected.csv
     */
    
    writeDetected("VarConjoin_detected.csv", stats, o);
    
    /*
     * Generating VarConjoin_linear.R
     */
    
    o.generate("VarConjoin_linear.R");
    o.writer->open("VarConjoin_linear.R");
    //o.writer->write(createVGROC("VarConjoin_detected.csv", "data$Depth", "'FP'"));
    o.writer->close();
}
