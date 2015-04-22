#include <iostream>
#include "standard.hpp"
#include "classify.hpp"
#include "abundance.hpp"
#include "expression.hpp"
#include "writers/r_writer.hpp"
#include "parsers/parser_tracking.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace SS::R;
using namespace Spike;

AbundanceStats Abundance::analyze(const std::string &file, const Options &options)
{
    AbundanceStats stats;
    const auto &r = Standard::instance();

    auto c = options.level == Gene ? countsForGenes() : countsForSequins();
    
    // Values for the x-axis and y-axis
    std::vector<double> x, y;

    ParserTracking::parse(file, [&](const Tracking &t)
    {
        assert(r.r_seqs_gA.count(t.geneID));

        switch (options.level)
        {
            case Gene:
            {
                c[t.geneID]++;
                const auto &m = r.r_seqs_gA.at(t.geneID);
                
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
                assert(r.r_seqs_iA.count(t.trackID));
                
                if (t.fpkm)
                {
                    const auto &i = r.r_seqs_iA.at(t.trackID);

                    x.push_back(i.abund(true));
                    y.push_back(t.fpkm);                    
                }
                
                break;
            }
        }
    });

    assert(!x.empty() && !y.empty());
    
    const auto cr = Expression::analyze(c);
    
    // Perform a linear-model to the abundance
    const auto m = lm(y, x); // correction bug

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

    if (options.level == Gene)
    {
        stats.sb = cr.sens(r.r_seqs_gA);
    }
    else
    {
        stats.sb = cr.sens(r.r_seqs_iA);
    }

    // Calculate the limit-of-sensitivity
    stats.sb = options.level == Gene ? cr.sens(r.r_seqs_gA) : cr.sens(r.r_seqs_iA);
    
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
                                                 % stats.sb.abund).str());
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