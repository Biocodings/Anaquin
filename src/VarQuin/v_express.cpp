#include "VarQuin/v_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotScatter();

VExpress::Stats VExpress::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    stats.n_gen = NAN;

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
            
            stats.n_syn++;
            stats.hist.at(d.id)++;            
        }
    });
    
    return stats;
}

void VExpress::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    /*
     * 1. Generating summary statistics
     */
    
    o.info("Generating VarExpress_summary.stats");
    o.writer->open("VarExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rAnnot,
                                                o.rAnnot,
                                                file,
                                                stats.hist,
                                                stats,
                                                stats,
                                                "sequins"));
    o.writer->close();
    
    /*
     * 2. Generating CSV for all sequins
     */
    
    o.info("Generating VarExpress_sequins.csv");
    o.writer->open("VarExpress_sequins.csv");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * 3. Generating for expression vs expression
     */
    
    o.info("Generating VarExpress_express.R");
    o.writer->open("VarExpress_express.R");
    o.writer->write(RWriter::createScript("VarExpress_sequins.csv", PlotScatter()));
    o.writer->close();
}