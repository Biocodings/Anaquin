#ifndef SS_TEST_HPP
#define SS_TEST_HPP

#include <ss/dist.hpp>
#include <ss/matrix.hpp>
#include <ss/internal/stest.hpp>

namespace SS
{
    template <typename T> static Results tTestOneSample
        (const T &x, Real h0, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::tTestOneSample(x, h0, conf, type);
    }
    
    /*
     * Two samples t-test, equal variance assumed.
     */

    template <typename Iter> static Results tTestTwoSamplesEqVar
            (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::tTestTwoSamples(x, y, h0, conf, type);
    }
    
    /*
     * Two samples paired t-test
     */

    template <typename Iter> static Results tTestPaired
            (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
    {
        return Internal::tTestPaired(x, y, h0, conf, type);
    }

    /*
     * Welch's t-test. An adaptation of Student's t-test and is more robust when the two
     * samples have unequal variances and unequal sample sizes (heteroscedastic).
     */

    template <typename Iter> static Results tTestWelch
            (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
    {
        return Internal::tTestWelch(x, y, h0, conf, type);
    }
}

#endif