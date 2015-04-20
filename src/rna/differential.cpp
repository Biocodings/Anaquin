#include "classify.hpp"
#include "expression.hpp"
#include "differential.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_cdiffs.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace SS::R;
using namespace Spike;

DifferentialStats Differential::analyze(const std::string &f, const Differential::Options &options)
{
    DifferentialStats stats;
    const auto &r = Standard::instance();

    auto c = options.level == Gene ? countsForGenes() : countsForSequins();

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
            case Gene:
            {
                assert(r.r_seqs_gA.count(t.geneID));
                assert(r.r_seqs_gB.count(t.geneID));

                if (t.status != NoTest)
                {
                    assert(t.fpkm_1);
                    assert(c.count(r.r_seqs_gA.at(t.geneID).geneID));
                    assert(c.count(r.r_seqs_gB.at(t.geneID).geneID));

                    // Calculate the known fold-change between B and A
                    known = r.r_seqs_gB.at(t.geneID).abund(true) / r.r_seqs_gA.at(t.geneID).abund(true);

                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;
                    
                    c[r.r_seqs_gA.at(t.geneID).geneID]++;
                    c[r.r_seqs_gB.at(t.geneID).geneID]++;
                    
                    x.push_back(known);
                    y.push_back(measured);
                }

                break;
            }

            case Isoform:
            {
                assert(r.r_seqs_iA.count(t.testID));
                assert(r.r_seqs_iB.count(t.testID));

                if (t.status != NoTest)
                {
                    // Calculate the known fold-change between B and A
                    known = r.r_seqs_iB.at(t.testID).raw / r.r_seqs_iA.at(t.testID).raw;

                    // Calculate the measured fold-change between B and A
                    measured = t.fpkm_2 / t.fpkm_1;

                    c[r.r_seqs_iA.at(t.testID).id]++;
                    c[r.r_seqs_iB.at(t.testID).id]++;

                    x.push_back(known);
                    y.push_back(measured);
                }

                break;
            }
        }
    });

    assert(!c.empty());
    assert(!x.empty() && x.size() == y.size());

    const auto cr = Expression::analyze(c);

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
    
    if (options.level == Gene)
    {
        stats.sb = cr.sens(r.r_seqs_gA); // TODO: Do we have to report for mixture B as well?
    }
    else
    {
        stats.sb = cr.sens(r.r_seqs_iA);
    }

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
    options.writer->write((boost::format(format) % stats.r2
                                                 % stats.slope
                                                 % stats.sb.abund).str());
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

    options.writer->write(RWriter::write(x, y));
    options.writer->close();

    return stats;
}