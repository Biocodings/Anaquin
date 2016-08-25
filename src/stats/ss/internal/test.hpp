#ifndef SS_INTERNAL_TEST_HPP
#define SS_INTERNAL_TEST_HPP

#include <functional>
#include <ss/stats.hpp>
#include <ss/data/errors.hpp>
#include <ss/data/results.hpp>
#include <ss/internal/dist.hpp>

namespace SS
{
    namespace Internal
    {
        inline Results logLRTest(Real l1, Real l2, DF d1, DF d2)
        {
            const auto lN = std::min(l1, l2);
            const auto lA = std::max(l1, l2);
            const auto dN = std::max(d1, d2);
            const auto dA = std::min(d1, d2);
            
            const auto df = dN - dA;
            const auto stats = -2.0 * (lN - lA);
            
            Results r;
            
            r.addDF(df);
            r.addStats(stats);
            r.addP(RMath::dchisq(stats, df, 0));
            
            return r;
        }
        
        inline Results chiSquareTest(const Matrix &exp, const Matrix &obs)
        {
            SS_ASSERT(exp.rows() == obs.rows() && exp.cols() == obs.cols(), "Dimenion mismatched");
            
            auto x2 = 0.0;
            
            MatrixUtils::each(exp, [&](Index i, Index j, Real ex)
            {
                const auto ob = obs(i, j);
                
                // Chi-square is sum of standard normals
                x2 += std::pow(ex - ob, 2) / ex;
            });
            
            const auto df = (exp.rows() - 1) * (exp.cols() - 1);
            const auto p  = 1.0 - RMath::pchisq(x2, df, 1, 0);
            
            Results r;
            
            r.addDF(df);
            r.addStats(x2);
            r.addP(p);
            
            return r;
        }

        inline Results chiSquareTestForIndependence(const Matrix &m)
        {
            // We'll need the marginal distribution
            const auto c = MatrixUtils::count(m);
            
            Matrix exp(m.rows(), m.cols());
            
            MatrixUtils::each(m, [&](Index i, Index j, Real observed)
            {
                const auto expected = (c.rsums(i) * c.csums(j) / c.sums);
                
                // The expected frequency under null hypothesis is the product of expected marginal distributions
                exp(i,j) = (c.rsums(i) * c.csums(j) / c.sums);
            });
            
            return chiSquareTest(exp, m);
        }
        
        inline Results mcNemarTest(const Matrix &m, bool yates = true)
        {
            SS_ASSERT(m.rows() == 2 && m.cols() == 2, "2x2 matrix required");
            
            const auto b = m(0,1);
            const auto c = m(1,0);
            const auto d = b-c-(yates ? 1.0 : 0);
            
            const auto df = 1;
            const auto x2 = (d*d) / (b+c);
            const auto p  = 1.0 - RMath::pchisq(x2, df, 1, 0);
            
            Results r;
            
            r.addStats(x2);
            r.addP(p);
            r.addDF(df);
            
            return r;
        }
        
        template <typename T> static Results chiSquareTestForVariance
                (const T &x, Real h0 = 1.0, P conf = 0.95, TestType type = TwoSided)
        {
            const auto n = x.size();
            const auto v = getVariance(x);
            
            const auto df = n - 1;
            const auto x2 = df * v / h0;
            
            const auto less    = RMath::pchisq(x2, df, 1, 0);
            const auto greater = 1.0 - less;
            
            auto alpha = conf.comp();
            
            if (type != TwoSided)
            {
                alpha = 2.0 * alpha;
            }
            
            const auto lc = type != Less ?    (df * v) / RMath::qchisq(1 - 0.5 * alpha, df, 1, 0) : 0.0;
            const auto uc = type != Greater ? (df * v) / RMath::qchisq(0.5 * alpha, df, 1, 0) : INFINITY;
            
            P p;
            
            switch (type)
            {
                case TestType::Less:     { p = less;    break; }
                case TestType::Greater:  { p = greater; break; }
                case TestType::TwoSided: { p = 2.0 * std::min(less, greater); break; }
            }
            
            Results r;
            
            r.addStats(x2);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
        
        template <typename Iter> static Results zTestOneSample
                (const Iter &x, Real h0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u = mean(x);
            const auto n = x.size();
            const auto s = sd(x) / sqrt(n);
            const auto z = (u - h0) / s;
            
            const auto c = Internal::critical(std::bind(Internal::qnorm, _1, 0, 1), conf, type);
            const auto p = Internal::pval(z, pnorm(z, 0, 1, 1, 0), type);
            
            const auto lc = (type != Less    ? u - (s * fabs(c)) : -INFINITY);
            const auto uc = (type != Greater ? u + (s * fabs(c)) :  INFINITY);
            
            Results r;
            
            r.addStats(z);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            
            return r;
        }
        
        template <typename Iter> static Results zTestTwoSamples
                (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            ASSERT(x.size() == y.size(), "Samples must match in dimension");
            ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u1 = mean(x);
            const auto u2 = mean(y);
            const auto s1 = sd(x);
            const auto s2 = sd(y);
            const auto n1 = x.size();
            const auto n2 = y.size();
            const auto u  = u1 - u2;
            
            const auto s = sqrt((((n1 - 1.0) * s1 * s1) + ((n2 - 1.0) * s2 * s2)) / (n1 + n2 - 2.0));
            
            const auto se = s * sqrt(1.0/n1 + 1.0/n2);
            
            const auto t = (u1 - u2 - h0) / se;
            const auto c = Internal::critical(std::bind(Internal::qnorm, _1, 0, 1), conf, type);
            const auto p = Internal::pval(t, pnorm(t, 0, 1, 1, 0), type);
            
            const auto lc = type != Less    ? (u1 - u2) - fabs(c) * se : -INFINITY;
            const auto uc = type != Greater ? (u1 - u2) + fabs(c) * se :  INFINITY;
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            
            return r;
        }
        
        /*
         * Two samples paired z-test
         */
        
        template <typename Iter> static Results zTestPaired
                (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            ASSERT(x.size() == y.size(), "Samples must match in dimension");
            ASSERT(x.size(), "Samples size must not be zero");
            
            std::vector<Real> diffs;
            diffs.resize(x.size());
            
            for (auto i = 0; i < x.size(); i++)
            {
                diffs[i] = x[i] - y[i];
            }
            
            return zTestOneSample(diffs, h0, conf, type);
        }
        
        /*
         * One sample z-test for proportion
         */
        
        inline Results zTestOneSampleProp
                (Counts x, Counts n, P p0 = 0.50, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(n, "Samples size must not be zero");
            SS_ASSERT(x <= n, "Number of success must be less than the sample size");
            
            const auto p = static_cast<Real>(x) / n;
            const auto s = sqrt(p0 * (1.0 - p0) / n);
            
            const auto z = (p - static_cast<Real>(p0)) / s;
            const auto c = Internal::critical(std::bind(Internal::qnorm, _1, 0, 1), conf, type);
            const auto v = Internal::pval(z, pnorm(z, 0, 1), type);
            
            const auto lc = (type != Less    ? p - (s * fabs(c)) : -INFINITY);
            const auto uc = (type != Greater ? p + (s * fabs(c)) :  INFINITY);
            
            Results r;
            
            r.addStats(z);
            r.addP(v);
            r.addLCI(lc);
            r.addUCI(uc);
            
            return r;
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
            using namespace std::placeholders;
            
            SS_ASSERT(n1 && n2, "Samples size must not be zero");
            SS_ASSERT(x1 <= n1 && x2 <= n2, "Number of success must be less than the sample size");
            
            const auto p1 = static_cast<Real>(x1) / n1;
            const auto p2 = static_cast<Real>(x2) / n2;
            
            // Estimate for the standard error
            const auto p = static_cast<Real>(x1 + x2) / (n1 + n2);
            
            // Uses the common estimate for the error
            const auto s = sqrt(p * (1.0 - p) * (1.0/n1 + 1.0/n2));
            
            const auto z = (p1 - p2 - static_cast<Real>(p0)) / s;
            const auto c = Internal::critical(std::bind(Internal::qnorm, _1, 0, 1), conf, type);
            const auto v = Internal::pval(z, pnorm(z, 0, 1), type);
            
            const auto lc = (type != Less    ? p1 - p2 - (s * fabs(c)) : -INFINITY);
            const auto uc = (type != Greater ? p1 - p2 + (s * fabs(c)) :  INFINITY);
            
            Results r;
            
            r.addStats(z);
            r.addP(v);
            r.addLCI(lc);
            r.addUCI(uc);
            
            return r;
        }

        template <typename Iter> static Results tTestOneSample
                (const Iter &x, Real h0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u  = getMean(x);
            const auto n  = x.size();
            const auto s  = getSD(x) / sqrt(n);
            const auto t  = (u - h0) / s;
            const auto df = n - 1;
            
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = (type != Less    ? u - (s * fabs(c)) : -INFINITY);
            const auto uc = (type != Greater ? u + (s * fabs(c)) :  INFINITY);
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
        
        template <typename Iter> static Results tTestTwoSamples
                (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u1 = getMean(x);
            const auto u2 = getMean(y);
            const auto s1 = getSD(x);
            const auto s2 = getSD(y);
            const auto n1 = x.size();
            const auto n2 = y.size();
            const auto u  = u1 - u2;
            
            const DF df = 2.0 * n1 - 2.0;
            
            const auto s  = sqrt((((n1 - 1.0) * s1 * s1) + ((n2 - 1.0) * s2 * s2)) / (n1 + n2 - 2.0));
            const auto se = s * sqrt(1.0/n1 + 1.0/n2);
            
            const auto t = (u1 - u2 - h0) / se;
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = type != Less    ? (u1 - u2) - fabs(c) * se : -INFINITY;
            const auto uc = type != Greater ? (u1 - u2) + fabs(c) * se :  INFINITY;
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
        
        template <typename Iter> static Results tTestPaired
        (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
        {
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            std::vector<Real> diffs;
            diffs.resize(x.size());
            
            for (auto i = 0; i < x.size(); i++)
            {
                diffs[i] = x[i] - y[i];
            }
            
            return tTestOneSample(diffs, h0, conf, type);
        }
        
        template <typename Iter> static Results tTestWelch
        (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u1 = getMean(x);
            const auto u2 = getMean(y);
            const auto s1 = getSD(x);
            const auto s2 = getSD(y);
            const auto n1 = x.size();
            const auto n2 = y.size();
            
            const auto s = std::sqrt((s1 * s1 / n1) + (s2 * s2 / n2));
            const auto t = (u1 - u2 - h0) / s;
            
            const DF df = std::pow(s1 * s1 / n1 + s2 * s2 / n2, 2) /
            (std::pow(s1 * s1 / n1, 2) / (n1 - 1) +
             std::pow(s2 * s2 / n2, 2) / (n2 - 1));
            
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = type != Less    ? (u1 - u2) - fabs(c) * s : -INFINITY;
            const auto uc = type != Greater ? (u1 - u2) + fabs(c) * s :  INFINITY;
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
    }
}

#endif