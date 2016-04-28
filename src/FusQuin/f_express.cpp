#include "FusQuin/f_express.hpp"
#include "parsers/parser_stab.hpp"
#include "parsers/parser_kallisto.hpp"

using namespace Anaquin;

// Defined resources.cpp
extern Scripts PlotFExpress();

FExpress::Stats FExpress::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_fus;
    
    FExpress::Stats stats;
    //stats.hist = r.fusionHist();

/*
    FUSQuin::analyze<FExpress::Options>(file, o, [&](const FUSQuin::Match &match)
    {
        switch (match.label)
        {
            case FUSQuin::Label::Positive:
            {
                stats.n_chrT++;
                break;
            }

            case FUSQuin::Label::GenoChrT:
            case FUSQuin::Label::Negative: { break; }

            case FUSQuin::Label::Geno:
            {
                stats.n_geno++;
                break;
            }
        }
    });
*/
    
    
    
    
    switch (o.caller)
    {
        case FusionCaller::StarFusion:
        {
            ParserSTab::parse(Reader(file), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
            {
                const SequinData *match;
                
                if ((match = r.findSplice(c.l)))
                {
                    const auto expected = match->mixes.at(Mix_1);
                    const auto measured = c.unique;

                    stats.add(match->id, expected, measured);
                    //stats.hist.at(match->id)++;
                }
            });

            break;
        }

        case FusionCaller::TopHatFusion:
        {
            throw "TopHatFusion not supported in FusExpress";
        }
    }

    stats.absolute = r.absolute(stats.hist);

    return stats;
}

void FExpress::report(const FileName &file, const Options &o)
{
    const auto stats = FExpress::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("FusExpress_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rGeno, file, stats.hist, stats, stats, ""));
    o.writer->close();

    /*
     * Generating CSV for all fusions
     */
    
    o.info("Generating FusExpress_quins.stats");
    o.writer->open("FusExpress_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating expression plot
     */
    
    o.writer->open("FusExpress_express.R");
    o.writer->write(RWriter::createScript("FusExpress_quins.stats", PlotFExpress()));
    o.writer->close();
}