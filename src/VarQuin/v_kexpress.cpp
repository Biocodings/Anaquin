#include "stats/kallisto.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_kexpress.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

// Defined by Kallisto
extern int __main__(int argc, char *argv[]);

VKExpress::Stats VKExpress::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    stats.n_endo = NAN;

    // Run quantification in Kallisto
    const auto abundFile = Kallisto::quant(o.index, file1, file2);
    
    /*
     * Parsing the generated files. We're interested in the file listing the abundance.
     */
    
    ParserKallisto::parse(Reader(abundFile), [&](const ParserKallisto::Data &d, const ParserProgress &)
    {
        const auto m = r.match(d.id);
        
        if (m)
        {
            // Expected abundance
            const auto known = m->mixes.at(Mix_1);
            
            // Measured abundance
            const auto measured = d.abund;
            
            stats.add(d.id, known, measured);
            
            stats.n_chrT++;
            stats.hist.at(d.id)++;            
        }
    });
    
    return stats;
}

void VKExpress::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating VarKExpress_summary.stats");
    o.writer->open("VarKExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                (file1 + " & " + file2),
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();
    
    /*
     * Generating CSV for all sequins
     */
    
    o.info("Generating VarKExpress_quins.csv");
    o.writer->open("VarKExpress_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();

    /*
     * Generating for AbundAbund
     */

    o.info("Generating VarKExpress_abundAbund.R");
    o.writer->open("VarKExpress_abundAbund.R");
    o.writer->write(RWriter::createScript("VarKExpress_quins.csv", PlotVAbundAbund()));
    o.writer->close();
}