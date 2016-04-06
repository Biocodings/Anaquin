#include "writers/r_writer.hpp"
#include "TransQuin/t_express.hpp"
#include "TransQuin/t_kexpress.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotTAbundAbund();

// Defined in resources.cpp
extern Scripts PlotRAbundAbund();

// Defined in resources.cpp
extern Scripts PlotMajor();

TKExpress::Stats TKExpress::analyze(const FileName &file, const Options &o)
{
    o.info("Parsing: " + file);

    const auto &r = Standard::instance().r_trans;
    
    TKExpress::Stats stats;
    
    ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &x, const ParserProgress &)
    {
        const TransData *m = nullptr;
        
        if ((m = r.match(x.id)))
        {
            stats.add(x.id, m->abund(Mix_1), x.abund);
        }
    });
    
    return stats;
}

void TKExpress::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TKExpress::analyze(files, o);
    
    o.info("Generating statistics");
    
    const auto units = "isoforms";
    
    /*
     * 1. Generating summary statistics (single or multiple samples)
     */
    
    o.info("Generating TransKExpress_summary.stats");
    o.writer->open("TransKExpress_summary.stats");
    
    if (files.size() == 1)
    {
        o.writer->write(TExpress::singleSummary(stats[0], files[0], units, o));
    }
    else
    {
        o.writer->write(TExpress::multipleSummary(files, stats, units, o));
    }

    o.writer->close();
    
    /*
     * 2. Generating detailed statistics for the sequins
     */
    
    o.info("Generating TransKExpress_quins.csv");
    o.writer->open("TransKExpress_quins.csv");
    
    if (files.size() == 1)
    {
        o.writer->write(StatsWriter::writeCSV(stats[0], "EAbund", "MAbund"));
    }
    else
    {
        o.writer->write(TExpress::multipleCSV(stats));
    }
    
    o.writer->close();
    
    /*
     * 3. Generating abundance vs abundance (single or multiple samples)
     */
    
    o.info("Generating TransKExpress_abundAbund.R");
    o.writer->open("TransKExpress_abundAbund.R");
    
    if (files.size() == 1)
    {
        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotTAbundAbund()));
    }
    else
    {
        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotRAbundAbund()));
    }
    
    o.writer->close();
    
    /*
     * 4. Generating major plot (but only if we have the isoforms...)
     */

    if (files.size() >= 2)
    {
        o.writer->open("TransKExpress_major.R");
        o.writer->write(RWriter::createScript("TransKExpress_quins.csv", PlotMajor()));
        o.writer->close();
    }
}