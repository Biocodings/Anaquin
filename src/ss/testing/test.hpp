#ifndef SS_TEST_HPP
#define SS_TEST_HPP

#include <ss/prob.hpp>
#include <ss/stats.hpp>

namespace SS
{
    enum TestType
    {
        OneTailed,
        TwoTailed,
    };

	template<typename T> struct Interval
	{
		T lower, upper;
	};

	template<typename T> struct TestResults
	{
        TestResults(TestStats x, const P p) : x(x), p(p) {}
        TestResults(TestStats x, DF df, const P p) : x(x), df(df), p(p) {}

        // Not all test has a degrees of freedom
        DF df = NAN;

        TestStats x;
		P p;
	};

	template <typename Dist> P statsTest(TestStats x, Dist dist, TestType type)
	{
		// Lowest significance level for which a null hypothesis still can be rejected
		const double cdf = P(dist.cdf(x));

		if (type == TwoTailed)
		{
			// The distribution is assumed be symmetry
			return P(2 * (cdf > 0.5 ? 1 - cdf : cdf));
		}
        else
        {
            // The distribution is assumed be symmetry
            return P(cdf > 0.5 ? 1 - cdf : cdf);
        }
	}
}

#endif