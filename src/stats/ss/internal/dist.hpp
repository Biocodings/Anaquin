#ifndef SS_INTERNAL_DIST_HPP
#define SS_INTERNAL_DIST_HPP

#include <math.h>
#include <ss/data/data.hpp>
#include <ss/data/errors.hpp>
#include <ss/internal/rmath.hpp>

namespace SS
{
    namespace Internal
    {
        inline double pt(double x, double n)
        {
            return RMath::pt(x, n);
        }
        
        inline double qt(double p, double ndf)
        {
            return RMath::qt(p, ndf, 1, 0);
        }
        
        inline double rnorm(double mu, double sigma)
        {
            return RMath::rnorm(mu, sigma);
        }
        
        inline double qnorm(double p, double mu, double sigma)
        {
            return RMath::qnorm(p, mu, sigma);
        }
        
        inline double dnorm(double x, double mu, double sigma)
        {
            return RMath::dnorm(x, mu, sigma);
        }
        
        inline double pnorm(double x, double mu, double sigma)
        {
            return RMath::pnorm(x, mu, sigma);
        }        
    }
}

#endif