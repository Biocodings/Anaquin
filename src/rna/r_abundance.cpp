#include <iostream>
#include "standard.hpp"
#include "classify.hpp"
#include "expression.hpp"
#include "r_abundance.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tracking.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace SS::R;
using namespace Spike;

RAbundanceStats RAbundance::analyze(const std::string &file, const Options &options)
{
    RAbundanceStats stats;
    const auto &s = Standard::instance();

    auto c = RAnalyzer::counter(options.level, options.mix);

    // Values for the x-axis and y-axis
    std::vector<double> x, y;

    ParserTracking::parse(file, [&](const Tracking &t, const ParserProgress &)
    {
        assert(s.r_seqs_gA.count(t.geneID));

        switch (options.level)
        {
            case Gene:
            {
                c[t.geneID]++;
                const auto &m = s.r_seqs_gA.at(t.geneID);
                
                if (t.fpkm)
                {
                    /*
                     * The x-axis would be the known concentration for each gene,
                     * the y-axis would be the expression (RPKM) reported.
                     */
                    
                    x.push_back(m.abund(true));
                    y.push_back(t.fpkm);
                }

                break;
            }

            case Isoform:
            {
                c[t.trackID]++;
                assert(s.r_seqs_iA.count(t.trackID));
                
                if (t.fpkm)
                {
                    const auto &i = s.r_seqs_iA.at(t.trackID);

                    x.push_back(i.abund(true));
                    y.push_back(t.fpkm);                    
                }
                
                break;
            }
        }
    });

    assert(!x.empty() && !y.empty());

    if (options.level == Gene)
    {
        stats.s = Expression::analyze(c, s.r_pair(options.mix));
    }
    else
    {
        stats.s = Expression::analyze(c, s.r_sequin(options.mix));
    }
    
    // Perform a linear-model to the abundance
    const auto m = lm(y, x);

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */

    stats.r2 = m.ar2;
    
    // Dependency between the two variables
    stats.r = cor(x, y);
    
    // Linear relationship between the two variables
    stats.slope = m.coeffs[1].value;

    const std::string format = "%1%\t%2%\t%3%";
    
    /*
     * Base-level statistics
     */
    
    if (options.level == Gene)
    {
        options.writer->open("abundance.genes.stats");
    }
    else
    {
        options.writer->open("abundance.isoform.stats");
    }

    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.r2
                                                 % stats.slope
                                                 % stats.s.abund).str());
    options.writer->close();

    /*
     * Generate a plot for the fold-change relationship
     */

    if (options.level == Gene)
    {
        options.writer->open("abundance.genes.R");
    }
    else
    {
        options.writer->open("abundance.isoform.R");
    }

    options.writer->write(RWriter::write(x, y));
    options.writer->close();

    return stats;
}