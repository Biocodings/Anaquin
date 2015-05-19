#ifndef SS_PROP_TEST_HPP
#define SS_PROP_TEST_HPP

#include <ss/dists/normal.hpp>
#include <ss/testing/test.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> struct ProportionTest
        {
            TestResults<T> test(unsigned n, unsigned x, T h) const
            {
                const auto p = P((double) x / n);
                const auto X = (p - h) / sqrt(h * (1 - h) / n);

                return TestResults<T>(X, statsTest(X, StandardNormal<T>(), TwoTailed));
            }
        };
    }

    namespace R
    {
        static Internal::ProportionTest<Real> prop;
    }
}

#endif