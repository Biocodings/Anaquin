#ifndef LINEAR_HPP
#define LINEAR_HPP

#include <map>
#include <vector>
#include "data/sample.hpp"

namespace Anaquin
{
    struct SLinearStats
    {
        SStrings files;
        
        // Pearson's correlation
        SReals r;
        
        // Regression slope
        SReals sl;
        
        // Coefficient of determination
        SReals R2;
        
        // F-statistic
        SReals F;
        
        // P-value under the null hypothesis
        SReals p;
        
        SReals SSM, SSE, SST;
        
        SCounts SSM_D, SSE_D, SST_D;
    };
    
    /*
     * Represents a simple linear regression fitted by maximum-likehihood estimation.
     *
     *   Model: y ~ c + m*x
     */
    
    struct LinearModel
    {
        // Constant coefficient
        double c;
        
        // Least-squared slope coefficient
        double m;
        
        // Adjusted R2
        double R2;
        
        // Pearson correlation
        double r;
        
        // Adjusted R2
        double aR2;
        
        double F, p;
        
        // ANOVA estimates
        double SST, SSM, SSE;
        
        // Degree of freedoms
        unsigned SST_D, SSM_D, SSE_D;
    };

    struct Point
    {
        Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        double x, y;
    };

    struct SequinStats : public std::map<SequinID, Point>
    {
        struct Data
        {
            std::vector<SequinID> ids;
            std::vector<double> x, y;

            std::map<SequinID, double> id2x;
            std::map<SequinID, double> id2y;
        };
        
        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }
        
        inline void sum(const SequinID &id, double x, double y)
        {
            assert((*this)[id].x == x);            
            (*this)[id].y += y;
        }
        
        Limit limitQuant() const;
        
        // Return the x-values and y-values after filtering
        Data data(bool shouldLog) const;
        
        // Compute a simple linear regression model. By default, this function assumes log-transformation.
        LinearModel linear(bool shouldLog = true) const;
    };
}

#endif
