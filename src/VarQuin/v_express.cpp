#include "VarQuin/v_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

VExpress::Stats VExpress::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    stats.n_endo = NAN;

    ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &d, const ParserProgress &)
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

void VExpress::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("VarExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT,
                                                o.rEndo,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();
    
    /*
     * Generating CSV for all sequins
     */
    
    o.writer->open("VarExpress_quins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating for AbundAbund
     */
    
    o.writer->open("VarExpress_abundAbund.R");
    o.writer->write(RWriter::createScript("VarExpress_quins.csv", PlotVAbundAbund()));
    o.writer->close();
}