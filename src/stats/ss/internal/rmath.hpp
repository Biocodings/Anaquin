#ifndef SS_INTERNAL4_HPP
#define SS_INTERNAL4_HPP

#include <vector>
#include <numeric>
#include <algorithm>
#include <ss/data/data.hpp>

namespace SS
{
    namespace RMath
    {
        extern "C" double rnorm(double mu, double sigma);
        extern "C" double dnorm4(double x, double mu, double sigma, int give_log);
        extern "C" double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);
        extern "C" double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);
        
        extern "C" double pf(double x, double df1, double df2, int lower_tail, int log_p);
        
        extern "C" double pt(double x, double n, int lower_tail, int log_p);
        extern "C" double qt(double p, double ndf, int lower_tail, int log_p);
        
        extern "C" double dchisq(double x, double df, int give_log);
        extern "C" double pchisq(double x, double df, int lower_tail, int log_p);
        extern "C" double qchisq(double p, double df, int lower_tail, int log_p);
        
        extern "C" double dsignrank(double x, double n, int give_log);
        extern "C" double psignrank(double x, double n, int lower_tail, int log_p);
        extern "C" double qsignrank(double x, double n, int lower_tail, int log_p);

        inline double pt(double x, double n) { return pt(x, n, 1, 0); }
        inline double qt(double x, double n) { return qt(x, n, 1, 0); }

        inline double rnorm(double x, double mu, double sigma) { return rnorm(mu, sigma); }
        inline double qnorm(double x, double mu, double sigma) { return qnorm5(x, mu, sigma, 1, 0); }
        inline double pnorm(double x, double mu, double sigma) { return pnorm5(x, mu, sigma, 1, 0); }
        inline double dnorm(double x, double mu, double sigma) { return dnorm4(x, mu, sigma, 0);    }
    }
}

#endif