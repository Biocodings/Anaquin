#include <iostream>
#include "standard.hpp"
#include "classify.hpp"
#include "abundance.hpp"
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
        assert(r.seqs_gA.count(t.geneID));
     
        switch (options.mode)
        {
            case AbdunanceGene:
            {
                c[t.geneID]++;
                const auto &m = r.seqs_gA.at(t.geneID);
                
                if (t.fpkm)
                {
                    /*
                     * The x-axis would be the known concentration for each gene,
                     * the y-axis would be the expression (RPKM) reported.
                     */
                    
                    x.push_back(m.r.exp + m.v.exp);
                    y.push_back(t.fpkm);
                }
                
                break;
            }

            case AbdunanceIsoform:
            {
                c[t.trackID]++;
                assert(r.seqs_iA.count(t.trackID));
                
                if (t.fpkm)
                {
                    const auto &i = r.seqs_iA.at(t.trackID);
                    
                    x.push_back(i.exp);
                    y.push_back(t.fpkm);                    
                }
                
                break;
            }
        }
    });

    const auto er = Expression::analyze(c);

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

    // Calcualte the limit-of-sensitivity
    stats.s = options.mode == AbdunanceGene ?
                Sensitivity(er.limit_key, er.limit_count,
                                r.seqs_gA.at(er.limit_key).r.exp + r.seqs_gA.at(er.limit_key).v.exp) :
                Sensitivity(er.limit_key, er.limit_count, r.seqs_iA.at(er.limit_key).exp);

    const std::string format = "%1%\t%2%\t%3%";
    
    options.writer->open("base.stats");
    options.writer->write((boost::format(format) % "r" % "s" % "ss").str());
    options.writer->write((boost::format(format) % stats.r2
                                                 % stats.slope
                                                 % stats.s.exp).str());
    options.writer->close();

    options.writer->open("scripts.R");
    options.writer->write(RWriter::write(y, x));
    options.writer->close();
    
    return stats;
}