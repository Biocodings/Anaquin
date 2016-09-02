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
    
    struct SInflectStats
    {
        SStrings files;

        // Eg: 2387648 (15.56%)
        SCounts countSyn, countGen;
        
        // Eg: sequins
        Units units;
        
        // Eg: Detected 46 sequins
        //SCounts n_det;
        
        // Optimal break-point
        SReals b;
        
        // Sequin for the break-point
        SStrings bID;
        
        // Regression intercepts
        SReals lInt, rInt;

        // Regression slopes
        SReals lSl, rSl;
        
        // Regression R2
        SReals lR2, rR2;
        
        // Pearson's correlation
        SReals lr, rr;
        
        // Linear regression with logarithm
        SLinearStats wLog;
    };

    struct LOQModel
    {
        // Name of the sequin
        SequinID id;
        
        // Pearson's correlation
        double lr = NAN, rr = NAN;
        
        // Coefficient of determination before and after the break-point
        double lR2 = NAN, rR2 = NAN;
        
        // Slope before and after the break-point
        double lSl = NAN, rSl = NAN;
        
        // Intercept before and after the break-point
        double lInt = NAN, rInt = NAN;

        // The optimal breakpoint
        double b = NAN;
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

    struct LinearStats : public std::map<SequinID, Point>
    {
        struct Data
        {
            std::vector<SequinID> ids;
            std::vector<double> x, y;

            std::map<SequinID, double> id2x;
            std::map<SequinID, double> id2y;
        };

        inline bool contains(const SequinID &id) const
        {
            return (*this).count(id);
        }
        
        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }
        
        inline void sum(const SequinID &id, double x, double y)
        {
            assert((*this)[id].x == x);            
            (*this)[id].y += y;
        }
        
        // Return the x-values and y-values after filtering
        Data data(bool shouldLog) const;
        
        // Compute the limit of quantification. By default, this function assumes log-transformation.
//        LOQModel limitQuant(bool shouldLog = true) const;
        
        // Compute a simple linear regression model. By default, this function assumes log-transformation.
        LinearModel linear(bool shouldLog = true) const;
    };
}

#endif
