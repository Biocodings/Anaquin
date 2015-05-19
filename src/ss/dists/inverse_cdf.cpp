#include <cmath>
#include <cassert>
#include <iostream>
#include "QStatistics/Distribution/InverseCDF.hpp"

using namespace QQ;

// Coefficients for the rational approximation.
const double InverseCumulativeNormal::a1_ = -3.969683028665376e+01;
const double InverseCumulativeNormal::a2_ =  2.209460984245205e+02;
const double InverseCumulativeNormal::a3_ = -2.759285104469687e+02;
const double InverseCumulativeNormal::a4_ =  1.383577518672690e+02;
const double InverseCumulativeNormal::a5_ = -3.066479806614716e+01;
const double InverseCumulativeNormal::a6_ =  2.506628277459239e+00;

const double InverseCumulativeNormal::b1_ = -5.447609879822406e+01;
const double InverseCumulativeNormal::b2_ =  1.615858368580409e+02;
const double InverseCumulativeNormal::b3_ = -1.556989798598866e+02;
const double InverseCumulativeNormal::b4_ =  6.680131188771972e+01;
const double InverseCumulativeNormal::b5_ = -1.328068155288572e+01;

const double InverseCumulativeNormal::c1_ = -7.784894002430293e-03;
const double InverseCumulativeNormal::c2_ = -3.223964580411365e-01;
const double InverseCumulativeNormal::c3_ = -2.400758277161838e+00;
const double InverseCumulativeNormal::c4_ = -2.549732539343734e+00;
const double InverseCumulativeNormal::c5_ =  4.374664141464968e+00;
const double InverseCumulativeNormal::c6_ =  2.938163982698783e+00;

const double InverseCumulativeNormal::d1_ =  7.784695709041462e-03;
const double InverseCumulativeNormal::d2_ =  3.224671290700398e-01;
const double InverseCumulativeNormal::d3_ =  2.445134137142996e+00;
const double InverseCumulativeNormal::d4_ =  3.754408661907416e+00;

#define QL_EPSILON             ((std::numeric_limits<double>::epsilon)())
#define QL_MIN_REAL           -((std::numeric_limits<double>::max)())
#define QL_MAX_REAL            ((std::numeric_limits<double>::max)())

const double InverseCumulativeNormal::x_low_ = 0.02425;
const double InverseCumulativeNormal::x_high_= 1.0 - x_low_;

    inline bool close_enough(double x, double y, std::size_t n) {
        // Deals with +infinity and -infinity representations etc.
        if (x == y)
            return true;

        double diff = std::fabs(x-y), tolerance = n * QL_EPSILON;

        if (x * y == 0.0) // x or y = 0.0
            return diff < (tolerance * tolerance);

        return diff <= tolerance*std::fabs(x) ||
               diff <= tolerance*std::fabs(y);
    }

	    inline bool close_enough(double x, double y) {
        return close_enough(x,y,42);
    }

double InverseCumulativeNormal::tail_value(double x) {
    if (x <= 0.0 || x >= 1.0) {
        // try to recover if due to numerical error
        if (close_enough(x, 1.0)) {
            return QL_MAX_REAL; // largest value available
        } else if (std::fabs(x) < QL_EPSILON) {
            return QL_MIN_REAL; // largest negative value available
        } else {
			assert(false);
            //QL_FAIL("InverseCumulativeNormal(" << x
            //        << ") undefined: must be 0 < x < 1");
        }
    }

    double z;
    if (x < x_low_) {
        // Rational approximation for the lower region 0<x<u_low
        z = std::sqrt(-2.0*std::log(x));
        z = (((((c1_*z+c2_)*z+c3_)*z+c4_)*z+c5_)*z+c6_) /
            ((((d1_*z+d2_)*z+d3_)*z+d4_)*z+1.0);
    } else {
        // Rational approximation for the upper region u_high<x<1
        z = std::sqrt(-2.0*std::log(1.0-x));
        z = -(((((c1_*z+c2_)*z+c3_)*z+c4_)*z+c5_)*z+c6_) /
            ((((d1_*z+d2_)*z+d3_)*z+d4_)*z+1.0);
    }

    return z;
}
