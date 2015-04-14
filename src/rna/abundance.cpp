#include <iostream>
#include "classify.hpp"
#include "abundance.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_tracking.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace Spike;

AbundanceStats Abundance::analyze(const std::string &file, const Abundance::AbundanceOptions &options)
{
    AbundanceStats stats;
    const auto r = StandardFactory::reference();

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
                
                /*
                 * The x-axis would be the known concentration for each gene,
                 * the y-axis would be the expression (RPKM) reported.
                 */
                
                x.push_back(m.r.exp + m.v.exp);
                y.push_back(t.fpkm);
                
                break;
            }

            case AbdunanceIsoform:
            {
                c[t.trackID]++;
                assert(r.seqs_iA.count(t.trackID));

                const auto &i = r.seqs_iA.at(t.trackID);

                x.push_back(i.exp);
                y.push_back(t.fpkm);
                
                break;
            }
        }
    });
    
    ANALYZE_COUNTS(c, cr);
    
    const auto lm = linearModel(y, x);

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */
    
    stats.r2 = lm.ar2;
    
    // Dependency between the two variables
    stats.r = pearson(x, y);
    
    // Linear relationship between the two variables
    stats.slope = lm.coeffs[1].value;

   // options.writer->write((boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%")
     //                      % "diluation" % "sn" % "sp" % "sensitivity").str());

    
    std::cout << stats.r2 << " " << stats.r << " " << stats.slope << std::endl;
    
    
    
    return stats;
}
