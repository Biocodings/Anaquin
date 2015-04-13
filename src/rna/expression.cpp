#include <limits>
#include <iostream>
#include "expression.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_ctracking.hpp"
#include <ss/regression/linear_model.hpp>

using namespace SS;
using namespace Spike;

ExpressionStats Expression::analyze(const std::string &file, ExpressionMode mode)
{
    const auto r = StandardFactory::reference();

    // Values for the x-axis and y-axis
    std::vector<double> x, y;

//    ParserCTracking::parse(file, [&](const CTracking &t)
//    {
//        assert(r.known(t.geneID));
//        assert(r.mix_gA.count(t.geneID));
//
//        switch (mode)
//        {
//            case GeneExpress:
//            {
//                const auto &a = r.mix_gA.at(t.geneID);
//
//                /*
//                 * The x-axis would be the known concentration for each gene, the y-axis would be the expression
//                 * (RPKM) reported.
//                 */
//
//                // The concentration for the gene is the sum of each isoform
//                //x.push_back(a.r.exp + a.v.exp);
//                
//                // The y-value is whatever reported
//                y.push_back(t.fpkm);
//                
//                break;
//            }
//
//            case IsoformExpress:
//            {
//                assert(r.mix_iA.count(t.trackID));
//                const auto &i = r.mix_iA.at(t.trackID);
//                
//                // The x-value is our known concentration
//                //x.push_back(i.exp);
//                
//                // The y-value is whatever reported
//                y.push_back(t.fpkm);
//                
//                break;
//            }
//        }
//    });
    
    const auto lm = linearModel(y, x);

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */
    
    ExpressionStats stats;
    
    stats.r2 = lm.ar2;
    
    // Dependency between the two variables
    stats.r = cor(x, y);
    
    // Linear relationship between the two variables
    stats.slope = lm.coeffs[1].value;

    std::cout << stats.r2 << " " << stats.r << " " << stats.slope << std::endl;
    
    return stats;
}
