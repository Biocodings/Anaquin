#include "FusQuin/f_fusion.hpp"
#include "parsers/parser_stab.hpp"

using namespace Anaquin;

// Defined resources.cpp
extern Scripts PlotFFusion();

FFusion::Stats FFusion::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_fus;
    
    FFusion::Stats stats;
    stats.hist = r.fusionHist();

    FUSQuin::analyze<FFusion::Options>(file, o, [&](const FUSQuin::Match &match)
    {
        switch (match.label)
        {
            case FUSQuin::Label::Positive:
            {
                const auto expected = r.match(match.known->id)->mixes.at(Mix_1);
                const auto measured = match.query.reads;
                
                stats.add(match.known->id, expected, measured);
                stats.hist.at(match.known->id)++;

                break;
            }
                
            case FUSQuin::Label::Geno:
            case FUSQuin::Label::GenoChrT:
            case FUSQuin::Label::Negative:
            {
                break;
            }
        }
    });
    
    stats.limit = r.absolute(stats.hist);

    return stats;
}

void FFusion::report(const FileName &file, const Options &o)
{
    const auto stats = FFusion::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("FusFusion_summary.stats");
    o.writer->write(StatsWriter::linearSummary(file, o.rAnnot, stats, stats, stats.hist, "sequins"));
    o.writer->close();

    /*
     * Generating CSV for all fusions
     */
    
    o.info("Generating FusFusion_sequins.stats");
    o.writer->open("FusFusion_sequins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating expression plot
     */
    
    o.info("Generating FusFusion_express.R");
    o.writer->open("FusFusion_express.R");
    o.writer->write(RWriter::createScript("FusFusion_sequins.stats", PlotFFusion()));
    o.writer->close();
}