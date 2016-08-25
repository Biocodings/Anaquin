#ifndef SS_TEST_HPP
#define SS_TEST_HPP

#include <ss/dist.hpp>
#include <ss/matrix.hpp>
#include <ss/internal/test.hpp>

namespace SS
{
    /*
     *                          Log Likelihood-Ratio Test
     *
     * Conduct the log likelihood-ratio test for two nested distribution models.
     *
     * When it is desired to incorporate covariates into an extreme value analysis, one method is to
     * incorporate them into the parameters of the extreme value distributions themselves in a regression-like
     * manner (cf. Coles, 2001 ch 6; Reiss and Thomas, 2007 ch.15). In order to justify whether or not
     * inclusion of the covariates into the model is significant or not is to apply the likelihood-ratio test
     * (of course, the test is more general than that, cf. Coles (2001) p.35).
     *
     * The test is only valid for comparing nested models. That is, the parameters of one model must be a
     * subset of the parameters of the second model.
     *
     * Suppose the base model, m0, is nested within the model m1. Let x be the negative log-likelihood for m0
     * and y for m1. Then the likelihood-ratio statistic (or deviance statistic) is given by (Coles, 2001,
     * p 35; Reiss and Thomas, 2007, p 118):
     *
     *         D = -2*(y - x).
     */
    
    /*
     * Performs LLRT for two hierarchical models.
     *
     *    - l1: Log-likelihood for the first model
     *    - l2: Log-likelihood for the second model
     *    - d1: Degree of freedom for the first model
     *    - d2: Degree of freedom for the second model
     */
    
    inline Results logLRTest(Real l1, Real l2, DF d1, DF d2)
    {
        return Internal::logLRTest(l1, l2, d1, d2);
    }
    
    /*
     *                      Pearson's Chi-Squared Test for Goodness of Fit
     *
     * The input matric are taken as two-dimensional contingency table: the entries must be non-negative
     * integers. The test is performed of the null hypothesis that the expected counts are consistent
     * to the measured counts.
     */
    
    inline Results chiSquareTest(const Matrix &exp, const Matrix &obs)
    {
        return Internal::chiSquareTest(exp, obs);
    }
    
    /*
     *                      Pearson's Chi-Squared Test For Independence
     *
     * The input matrix is taken as a two-dimensional contingency table: the entries must be non-negative
     * integers. The test is performed of the null hypothesis that the joint distribution of the cell counts
     * in a 2-dimensional contingency table is the product of the row and column marginals.
     */
    
    inline Results chiSquareTestForIndependence(const Matrix &m)
    {
        return Internal::chiSquareTestForIndependence(m);
    }
    
    /*
     * -------------------- McNemar's Test --------------------
     */
    
    inline Results mcNemarTest(const Matrix &m, bool yates = true)
    {
        return Internal::mcNemarTest(m, yates);
    }
    
    /*
     * -------------------- Pearson's Chi-Squared Test For Variance --------------------
     */
    
    template <typename T> static Results chiSquareTestForVariance
            (const T &x, Real h0 = 1.0, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::chiSquareTestForVariance(x, h0, conf, type);
    }
    
    template <typename T> static Results zTestOneSample
        (const T &x, Real h0, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::zTestOneSample(x, h0, conf, type);
    }
    
    template <typename T> static Results zTestTwoSamples
        (const T &x, const T &y, Real h0 = 0.0, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::zTestTwoSamples(x, y, h0, conf, type);
    }
    
    /*
     * Two samples paired z-test
     */
    
    template <typename T> static Results zTestPaired
        (const T &x, const T &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
    {
        return Internal::zTestPaired(x, y, h0, conf, type);
    }
    
    /*
     * One sample z-test for proportion
     */
    
    inline Results zTestOneSampleProp
        (Counts x, Counts n, P p0 = 0.50, P conf = 0.95, TestType type = TwoSided)
    {
        return Internal::zTestOneSampleProp(x, n, p0, conf, type);
    }
    
    /*
     * Two samples z-test for proportion
     */
    
    inline Results zTestTwoSamplesProp(Counts x1,
                                       Counts n1,
                                       Counts x2,
                                       Counts n2,
                                       P p0 = 0.00,
                                       P conf = 0.95,
                                       TestType type = TwoSided)
    {
        return Internal::zTestTwoSamplesProp(x1, n1, x2, n2, p0, conf, type);
    }
    
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