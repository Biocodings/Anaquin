#include "stats/kallisto.hpp"
#include "VarQuin/v_kallele.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotAlleleAllele();

VKAllele::Stats VKAllele::analyze(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VKAllele::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();
    
    stats.n_endo = NAN;
    
    // Run quantification in Kallisto
    Kallisto::quant(o.file, file1, file2);

    /*
     * Parse the generated files. We're interested in the file listing the abundance.
     */
    
    ParserKallisto::parse(Reader(Kallisto::abundFile), [&](const ParserKallisto::Data &d, const ParserProgress &)
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
    
    //stats.all.limit = r.absolute(stats.hist);
    
    return stats;
}

void VKAllele::report(const FileName &file1, const FileName &file2, const Options &o)
{
    const auto &stats = analyze(file1, file2, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */

    o.info("Generating VarKAllele_summary.stats");
    o.writer->open("VarKAllele_summary.stats");
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

    o.info("Generating VarKAllele_quins.csv");
    o.writer->open("VarKAllele_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating for AlleleAllele
     */

    o.info("Generating VarKAllele_alleleAllele.R");
    o.writer->open("VarKAllele_alleleAllele.R");
    o.writer->write(RWriter::createScript("VarKAllele_quins.csv", PlotAlleleAllele()));
    o.writer->close();
}