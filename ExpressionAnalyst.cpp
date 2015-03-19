#include "StandardFactory.hpp"
#include "ParserCTracking.hpp"
#include "ExpressionAnalyst.hpp"
#include "Stats/Regression/LinearRegression.hpp"

using namespace QQ;

ExpressionStats ExpressionAnalyst::analyze(const std::string &file, ExpressionMode mode, Sequins s, Reads n)
{
    auto r = StandardFactory::reference();

    // Values for the x-axis and y-axis
    std::vector<double> x, y;

    ParserCTracking::parse(file, [&](const CTracking &t)
    {
        assert(r.known(t.id));
        assert(r.mixA.count(t.id));

        switch (mode)
        {
            case GeneExpress:
            {
                const auto &a = r.mixA[t.id];

                /*
                 * The x-axis would be the known concentration for each gene, the y-axis would be the expression
                 * (RPKM) reported.
                 */

                // The concentration for the gene is the sum of each isoform
                x.push_back(a.r.exp + a.v.exp);
                
                // The y-value is whatever reported
                y.push_back(t.fpkm);
                
                break;
            }

            case IsoformsExpress:
            {
                assert(r.isoA.count(t.id));
                
                const auto &i = r.isoA[t.id];
                
                // The x-value is our known concentration
                x.push_back(i.exp);
                
                // The y-value is whatever reported
                y.push_back(t.fpkm);
                
                break;
            }
        }
    });
    
    const auto lm = linearModel(y, x);

    /*
     * In our analysis, the dependent variable is expression while the independent
     * variable is the known concentraion.
     *
     *     expression = constant + slope * concentraion
     */
    
    ExpressionStats stats;
    
    stats.r2 = lm.ar2;
    
    // Bounded correlation between the two variables
    stats.r = cor(x, y);
    
    // Estimated change in the expected value for to a 1-unit increase
    stats.slope = lm.coeffs[1].value;

    std::cout << stats.r2 << " " << stats.r << " " << stats.slope << std::endl;
    
	return stats;
}