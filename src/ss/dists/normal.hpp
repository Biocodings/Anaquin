#ifndef SS_NORMAL_HPP
#define SS_NORMAL_HPP

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include <ss/prob.hpp>
#include <ss/stats.hpp>
#include <ss/types.hpp>
#include <ss/dists/inverse_cdf.hpp>
#include <ss/dists/mersenne_twister.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> struct StandardNormal
        {
            P cdf(double x) const;

            // Generate a pseudo-random standard normal
            T random();
        };

        template <typename T> struct Normal : public StandardNormal<T>
        {
            Normal<T>(T mean, T sd);
            
            P cdf(double x) const;
            
            // Generate a pseudo-random normal
            T random();

            T mean, sd;
        };

        template <typename T> P StandardNormal<T>::cdf(double x) const
        {
            const T k = 1.0 / (1.0 + 0.2316419*x);
            const T k_sum = k * (0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

            if (x >= 0.0)
            {
                return P(1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
            }
            else
            {
                return P(1.0 - cdf(-x));
            }
        }

        template <typename T> Normal<T>::Normal(T mean, T sd) : mean(mean), sd(sd) {}

        template <typename T> T random()
        {
            static InverseCDF<MersenneTwisterUniformRng, InverseCumulativeNormal, Real> g;
            return g();
        }
        
        template <typename T> P Normal<T>::cdf(double x) const
        {
            return StandardNormal<T>::cdf((x - mean) / sd);
        }

        template <typename T> T Normal<T>::random()
        {
            return mean + StandardNormal<T>::random() * sd;
        }
    }

    typedef Internal::Normal<Real> Normal;
    typedef Internal::StandardNormal<Real> StandardNormal;

	inline P N(double x)
	{
		return StandardNormal().cdf(x);
	}
}

#endif