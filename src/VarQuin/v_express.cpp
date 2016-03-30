#include "VarQuin/v_express.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVAbundAbund();

static void writeCSV(const FileName &file, const VExpress::Stats &stats, const VExpress::Options &o)
{
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%";
    
    auto f = [&](const LinearStats &l)
    {
        const auto data = l.data(false);
        
        for (auto i = 0; i < data.ids.size(); i++)
        {
            o.writer->write((boost::format(format) % data.ids[i]
                                                   % data.x[i]
                                                   % data.y[i]).str());
        }
    };
    
    o.writer->write((boost::format(format) % "Sequin"
                                           % "EAbund"
                                           % "MAbund").str());
    f(stats);

    o.writer->close();
}

VExpress::Stats VExpress::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VExpress::Stats stats;
    
    // Initialize the distribution for each sequin
    stats.hist = r.hist();

    ParserKallisto::parse(Reader(file), [&](const ParserKallisto::Data &d, const ParserProgress &)
    {
        const auto m = r.match(d.id);
        
        if (m)
        {
            // Expected allele frequency
            const auto known = m->mixes.at(Mix_1); // alleleFreq(m->id);

            // Measured abundance
            const auto measured = d.abund;
            
            stats.add(d.id, known, measured);
        }
        else
        {
            std::cout << "d" << std::endl;
        }
    });
    
    return stats;
}

void VExpress::report(const FileName &file, const Options &o)
{
    const auto &stats = analyze(file, o);
    
    /*
     * Generating CSV for all sequins
     */
    
    writeCSV("VarExpress_quins.csv", stats, o);
    
    /*
     * Generating for AlleleAllele
     */
    
    o.writer->open("VarExpress_abundAbund.R");
    o.writer->write(RWriter::createScript("VarExpress_quins.csv", PlotVAbundAbund()));
    o.writer->close();
    
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
                                                "variants"));
    o.writer->close();
}