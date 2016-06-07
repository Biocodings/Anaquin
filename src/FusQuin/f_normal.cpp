#include "FusQuin/f_normal.hpp"
#include "parsers/parser_stab.hpp"

using namespace Anaquin;

// Defined resources.cpp
extern Scripts PlotFNormal();

FNormal::Stats FNormal::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_fus;
    
    FNormal::Stats stats;
    stats.hist = r.normalHist();

    switch (o.soft)
    {
        case FusionCaller::StarFusion:
        {
            ParserSTab::parse(Reader(file), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
            {
                const SequinData *match;
                
                if (c.id == ChrT) { stats.n_chrT++;
                }
                else
                {
                    stats.n_geno++;
                }
                
                if (c.id == ChrT && (match = r.findJunct(c.l)))
                {
                    const auto expected = match->mixes.at(Mix_1);
                    const auto measured = c.unique;
                    
                    stats.add(match->id, expected, measured);
                    stats.hist.at(match->id)++;
                }
            });

            break;
        }

        case FusionCaller::TopHatFusion:
        {
            throw "TopHatFusion not supported in FusNormal";
        }
    }

    stats.limit = r.absolute(stats.hist);

    return stats;
}

void FNormal::report(const FileName &file, const Options &o)
{
    const auto stats = FNormal::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating FusNormal_summary.stats
     */
    
    o.generate("FusNormal_summary.stats");
    o.writer->open("FusNormal_summary.stats");
    o.writer->write(StatsWriter::linearSummary(file, o.rAnnot, stats, stats, stats.hist, "sequins"));
    o.writer->close();

    /*
     * Generating FusNormal_quins.stats
     */
    
    o.generate("FusNormal_quins.stats");
    o.writer->open("FusNormal_quins.stats");
    o.writer->write(StatsWriter::writeCSV(stats));
    o.writer->close();
    
    /*
     * Generating FusNormal_express.R
     */
    
    o.generate("FusNormal_express.R");
    o.writer->open("FusNormal_express.R");
    o.writer->write(RWriter::createScript("FusNormal_quins.stats", PlotFNormal()));
    o.writer->close();
}