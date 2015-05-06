#include <ss/c.hpp>
#include "expression.hpp"
#include "r_differential.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace Spike;

RDifferentialStats RDifferential::analyze(const std::string &f, const Options &options)
{
    RDifferentialStats stats;
    const auto &s = Standard::instance();

    auto c = (options.level == Gene ? RAnalyzer::geneCounter() : RAnalyzer::isoformCounter());

    std::vector<Concentration> x, y;
    std::vector<std::string> z;

    ParserCDiffs::parse(f, [&](const TrackingDiffs &t, const ParserProgress &)
    {
        Fold known, measured;

        /*
         * By measuring and comparing the observed fold-changes with the known changes,
         * it's possible to create a linear model. In a perfect experiment, one would
         * expect perfect correlation.
         *
         * For example, let's say R_1_1 is a silico gene. We might have the following table
         *
         *        R_1_1, 10000000, 2500000
         *
         * This is a fold-change of 2.5. One would expect a similar fold-change observed
         * in the experiment.
         */

        switch (options.level)
        {
            case Gene:
            {
                if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    // Calculate the known fold-change between B and A
                    known = s.r_seqs_gB.at(t.geneID).abund(true) / s.r_seqs_gA.at(t.geneID).abund(true);

                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;

                    c[t.geneID]++;                    
                    x.push_back(known);
                    y.push_back(measured);
                    z.push_back(t.geneID);
                }

                break;
            }

            case Isoform:
            {
               if (t.status != NoTest && t.fpkm_1 && t.fpkm_2)
                {
                    // Calculate the known fold-change between B and A
                    known = s.r_seqs_iB.at(t.testID).raw / s.r_seqs_iA.at(t.testID).raw;

                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;

                    if (known)
                    {
                        c[t.testID]++;
                        x.push_back(known);
                        y.push_back(measured);
                        z.push_back(t.testID);
                    }
                }

                break;
            }
        }
    });

    assert(!c.empty() && !x.empty());
    assert(!x.empty() && x.size() == y.size());

    if (options.level == Gene)
    {
        stats.s = Expression::analyze(c, s.r_pair(options.rMix));
    }
    else
    {
        stats.s = Expression::analyze(c, s.r_sequin(options.rMix));
    }

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */

    const auto m = lm("y ~ x", data.frame(SS::c(y), SS::c(x)));

    stats.lm.r2 = m.ar2;
    stats.lm.r = cor(x, y);
    stats.lm.m = m.coeffs[1].v;

    /*
     * Base-level statistics
     */
    
    const std::string format = "%1%\t%2%\t%3%";

    if (options.level == Gene)
    {
        options.writer->open("diffs.genes.stats");
    }
    else
    {
        options.writer->open("diffs.isoform.stats");
    }
    
    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.lm.r
                                                 % stats.lm.m
                                                 % stats.s.abund).str());
    options.writer->close();

    /*
     * Generate a plot for the fold-change relationship
     */
    
    if (options.level == Gene)
    {
        options.writer->open("diffs.genes.R");
    }
    else
    {
        options.writer->open("diffs.isoform.R");
    }

    options.writer->write(RWriter::write(x, y, z));
    options.writer->close();

    return stats;
}