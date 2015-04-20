#include "classify.hpp"
#include "expression.hpp"
#include "differential.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace SS::R;
using namespace Spike;

DifferentialStats Differential::analyze(const std::string &f, const Differential::Options &options)
{
    DifferentialStats stats;
    const auto &r = Standard::instance();

    INIT_COUNTER(c);

    // Values for the coordinates
    std::vector<Concentration> x, y;

    ParserCDiffs::parse(f, [&](const TrackingDiffs &t)
    {
        Fold known, measured;

        /*
         * In a differential-expression experiment, mixtures of A and B would be spiked into two samples
         * respectively. By measuring and comparing the observed fold-changes with the known changes,
         * it's possible to perform a linear-regression model. In a perfect experiment, one would
         * expect perfect correlation.
         *
         * For example, let's say R_1_1 is a silico gene. We might have the following table
         *
         *        R_1_1, 10000000, 2500000
         *
         * This is a fold-change of 2.5. One would expect a similar fold-change observed
         * in the experiment. Please refer to the documentation for more details.
         */

        switch (options.level)
        {
            case LevelGene:
            {
                assert(r.seqs_gA.count(t.geneID));
                assert(r.seqs_gB.count(t.geneID));

                if (t.status != NoTest)
                {
                    assert(t.fpkm_1);
                    
                    // Calculate the known fold-change between B and A
                    known = r.seqs_gB.at(t.geneID).raw() / r.seqs_gA.at(t.geneID).raw();
                    
                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;
                    
                    c[r.seqs_gA.at(t.geneID).r.id]++;
                    c[r.seqs_gA.at(t.geneID).v.id]++;
                    c[r.seqs_gB.at(t.geneID).r.id]++;
                    c[r.seqs_gB.at(t.geneID).v.id]++;
                    
                    x.push_back(known);
                    y.push_back(measured);
                }

                break;
            }

            case LevelIsoform:
            {
                assert(r.seqs_iA.count(t.testID));
                assert(r.seqs_iB.count(t.testID));

                if (t.status != NoTest)
                {
                    // Calculate the known fold-change between B and A
                    known = r.seqs_iB.at(t.testID).raw / r.seqs_iA.at(t.testID).raw;

                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;

                    c[r.seqs_iA.at(t.testID).id]++;
                    c[r.seqs_iB.at(t.testID).id]++;

                    x.push_back(known);
                    y.push_back(measured);
                }

                break;
            }
        }
    });

    assert(!c.empty());
    assert(!x.empty() && x.size() == y.size());

    const auto er = Expression::analyze(c);

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */
    
    const auto m = lm(y, x);

    stats.r2 = m.ar2;
    
    // Dependency between the two variables
    stats.r = cor(x, y);

    // Linear relationship between the two variables
    stats.slope = m.coeffs[1].value;
    
    stats.s = er.sens();

    const std::string format = "%1%\t%2%\t%3%";

    options.writer->open("base.stats"); // Name???
    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.r2
                                                 % stats.slope
                                                 % stats.s.exp).str());
    options.writer->close();

    return stats;
}