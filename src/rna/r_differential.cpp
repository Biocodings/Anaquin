#include "expression.hpp"
#include "r_differential.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace Spike;

RDifferentialStats RDifferential::analyze(const std::string &f, const Options &options)
{
    RDifferentialStats stats;
    const auto &s = Standard::instance();

    auto c = (options.level == Gene ? RAnalyzer::geneCounter() : RAnalyzer::isoformCounter());

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
                    stats.x.push_back(log(known));
                    stats.y.push_back(log(measured));
                    stats.z.push_back(t.geneID);
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
                        stats.x.push_back(log(known));
                        stats.y.push_back(log(measured));
                        stats.z.push_back(t.testID);
                    }
                }

                break;
            }
        }
    });

    assert(!c.empty() && !stats.x.empty());
    assert(!stats.x.empty() && stats.x.size() == stats.y.size());

    stats.linear();

    if (options.level == Gene)
    {
        stats.s = Expression::analyze(c, s.r_pair(options.rMix));
        AnalyzeReporter::report("diffs.stats", "diffs.genes.R", stats, c, options.writer);
    }
    else
    {
        stats.s = Expression::analyze(c, s.r_sequin(options.rMix));
        AnalyzeReporter::report("diffs.stats", "diffs.isoform.R", stats, c, options.writer);
    }

    return stats;
}