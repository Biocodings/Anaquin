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

AbundanceStats Abundance::analyze(const std::string &file, const Abundance::Options &options)
{
    AbundanceStats stats;
    const auto &r = Standard::instance();

    INIT_COUNTER(c);
    
    // Values for the x-axis and y-axis
    std::vector<double> x, y;

    ParserTracking::parse(file, [&](const Tracking &t)
    {
        assert(r.r_seqs_gA.count(t.geneID));

        switch (options.level)
        {
            case LevelGene:
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

            case LevelIsoform:
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
    
    const auto er = Expression::analyze(c);
    
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

    /*
    // Calculate the limit-of-sensitivity
    stats.s = options.level == LevelGene ?
                  Sensitivity(er.limit_key, er.limit_count,
                                r.r_seqs_gA.at(er.limit_key).r.raw + r.r_seqs_gA.at(er.limit_key).v.raw) :
                  Sensitivity(er.limit_key, er.limit_count, r.r_seqs_iA.at(er.limit_key).raw);
     */
    
    const std::string format = "%1%\t%2%\t%3%";
    
    options.writer->open("base.stats"); // should be named as gene_exp.stat or isoform_exp.stats
    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.r2
                                                 % stats.slope
                                                 % stats.s.abund).str());
    options.writer->close();

    options.writer->open("scripts.R");
    options.writer->write(RWriter::write(y, x));
    options.writer->close();
    
    return stats;
}