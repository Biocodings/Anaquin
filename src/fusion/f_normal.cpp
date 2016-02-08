#include "fusion/f_normal.hpp"
#include "parsers/parser_stab.hpp"

using namespace Anaquin;

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

                    stats.data.add(match->id, expected, measured);
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

void FNormal::report(const FileName &file, const Options &o)
{
    const auto stats = FNormal::analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating summary statistics
     */
    
    o.writer->open("FusionNormal_summary.stats");
    o.writer->write(StatsWriter::inflectSummary(o.rChrT(), o.rEndo(), file, stats, stats.data, ""));
    o.writer->close();

    /*
     * Generating scatter plot
     */
    
    //AnalyzeReporter::scatter(stats, "", "FusionNormal", "Expected concentration (attomol/ul)", "Measured coverage (reads)", "Expected concentration (log2 attomol/ul)", "Measured coverage (log2 reads)", o.writer);
}