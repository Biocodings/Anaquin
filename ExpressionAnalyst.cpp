#include "StandardFactory.hpp"
#include "ParserCTracking.hpp"
#include "ExpressionAnalyst.hpp"
#include "Stats/Regression/LinearRegression.hpp"

using namespace QQ;

ExpressionStats ExpressionAnalyst::analyze(const std::string &file, Sequins s, Reads n)
{
    const auto r = StandardFactory::reference();

    // Values for the x-axis and y-axis
    std::vector<float> x, y;
    
    ParserCTracking::parse(file, [&](const CTracking &t)
    {
        assert(r.knownGene(t.geneID));
        const auto c = r.concent(t.geneID);
        
        x.push_back(c.amounts);
        y.push_back(t.fpkm);
    });
    
    const auto lm = linearModel(y, x);

    /*
     * In our analysis, the dependent variable is expression (FPKM) while the independent
     * variable is known concentraion.
     *
     *      expression = constant + slope * concentraion
     */
    
    ExpressionStats stats;
    
    stats.r2 = lm.ar2;
    
    // Bounded correlation between the two variables
    stats.r = cor(x, y);
    
    // Estimated change in the expected value for to a 1-unit increase
    stats.slope = lm.coeffs[1].value;

	return stats;
}