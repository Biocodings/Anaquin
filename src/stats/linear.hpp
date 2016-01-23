#ifndef LINEAR_HPP
#define LINEAR_HPP

#include <map>
#include <vector>
#include "stats/limit.hpp"
#include "stats/stats.hpp"

namespace Anaquin
{
    struct SLinearStats
    {
        SStrings files;
        
        // Pearson's correlation
        SReals cor;
        
        // Regression slope
        SReals slope;
        
        // Coefficient of determination
        SReals R2;
        
        // F-statistic
        SReals F;
        
        // P-value under the null hypothesis
        SReals p;
        
        SReals SSM, SSE, SST;
        
        SCounts SSM_D, SSE_D, SST_D;
    };
    
    struct SInflectStats
    {
        SStrings files;
        
        // Eg: 2387648 (15.56%)
        SReals chrT_n, endo_n;
        
        // Eg: 25738262 (67.56%)
        SProps chrT_p, endo_p;
        
        // Eg: 56 sequins
        SCounts ref_n;
        
        // Eg: sequins
        Units units;
        
        // Eg: Detected 46 sequins
        SCounts det_n;
        
        // Optimal break-point
        SReals b;
        
        // Sequins for the break-points
        SStrings bID;
        
        // Regression intercepts
        SReals lInter, rInter;
        
        // Regression slopes
        SReals lSlope, rSlope;
        
        // Regression R2
        SReals lR2, rR2;
        
        // Linear regression without logarithm
        SLinearStats nLog;
        
        // Linear regression with logarithm
        SLinearStats wLog;
    };
    

    
    
    
    struct InflectionLimit
    {
        // Name of the sequin
        std::string id;
        
        // Coefficient of determination before and after the break-point
        double lR2, rR2;
        
        // Slope before and after the break-point
        double lSl, rSl;
        
        // Intercept before and after the break-point
        double lInt, rInt;
        
        // The optimal breakpoint
        double b;
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
        double r2;
        
        // Pearson correlation
        double r;
        
        // Adjusted R2
        double ar2;
        
        double f, p;
        double sst, ssm, sse;
        
        // Degree of freedoms
        unsigned sst_df, ssm_df, sse_df;
    };

    struct Point
    {
        Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        double x, y;
    };

    struct LinearStats : public std::map<SequinID, Point>
    {
        Limit s;
        
        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }
        
        // Return the x-values and y-values after filtering
        void data(std::vector<double> &x, std::vector<double> &y, bool shouldLog, std::vector<std::string> *ids = nullptr) const;
        
        // Compute the inflection limit. By default, this function assumes log-transformation.
        InflectionLimit inflect(bool shouldLog = true) const;
        
        // Compute a simple linear regression model. By default, this function assumes log-transformation.
        LinearModel linear(bool shouldLog = true) const;
    };
}

#endif
