#include "fusion/f_normal.hpp"
#include "parsers/parser_stab.hpp"

using namespace Anaquin;

// Defined resources.cpp
extern Scripts PlotNormal();

FNormal::Stats FNormal::analyze(const FileName &splice, const Options &o)
{
    FNormal::Stats stats;
    
    const auto &r = Standard::instance().r_fus;

    switch (o.caller)
    {
        case FusionCaller::Star:
        {
            ParserSTab::parse(Reader(splice), [&](const ParserSTab::Chimeric &c, const ParserProgress &)
            {
                if (c.id == ChrT)
                {
                    stats.n_chrT++;
                }
                else
                {
                    stats.n_endo++;
                }

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

        case FusionCaller::TopHat:
        {
            throw "TopHat not supported in FusionNormal";
        }
    }

    stats.limit = r.limit(stats.hist);

    return stats;
}

static void writeCSV(const FileName &file, const FNormal::Stats &stats, const FNormal::Options &o)
{
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->write((boost::format(format) % "Sequin"
                                           % "EAbund"
                                           % "MAbund").str());
    
    const auto data = stats.data(false);

    for (auto i = 0; i < data.ids.size(); i++)
    {
        o.writer->write((boost::format(format) % data.ids[i]
                                               % data.x[i]
                                               % data.y[i]).str());
    }
    
    o.writer->close();
}

void FNormal::report(const FileName &file, const Options &o)
{
    const auto stats = FNormal::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("FusionNormal_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT, o.rEndo, file, stats.hist, stats, stats, ""));
    o.writer->close();

    /*
     * Generating CSV for all fusions
     */
    
    writeCSV("FusionNormal_quins.csv", stats, o);

    /*
     * Generating scatter plot
     */
    
    o.writer->open("FusionNormal_scatter.R");
    o.writer->write(RWriter::createScript("FusionNormal_quins.csv", PlotNormal()));
    o.writer->close();
}