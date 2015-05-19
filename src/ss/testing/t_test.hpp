#ifndef SS_T_TEST_HPP
#define SS_T_TEST_HPP

#include <ss/dists/t_dist.hpp>
#include <ss/testing/test.hpp>

namespace SS
{
    /*
     * Unpaired one-group t-test
     */
    
    template <typename T, typename Iter> static TestResults<T> tt_one_group
            (const Iter &x, TestType type = TwoTailed, T h0 = T())
    {
        const auto u  = mean(x);
        const auto s  = sd(x);
        const auto n  = x.size();
        const auto t  = (u - h0) / sqrt(s / n);
        const auto df = n - 1;

        return TestResults<T>(t, df, statsTest<TestResults<T>>(t, TDist<T>(df), type));
    }

    /*
     * Unpaired two-groups t-test
     */

    template <typename T, typename Iter> static TestResults<T> tt_two_groups
            (const Iter &x, const Iter &y, TestType type = TwoTailed, T h0 = T())
    {
        const auto u1 = mean(x);
        const auto u2 = mean(y);
        const auto s1 = sd(x);
        const auto s2 = sd(y);
        const auto n1 = x.size();
        const auto n2 = y.size();
        const auto s  = sqrt((((n1 - 1) * s1 * s1) + ((n2 - 1) * s2 * s2)) / (n1 + n2 - 2));
        const auto t  = (u1 - u2 - h0) / (s * sqrt(1.0/n1 + 1.0/n2));
        const DF df   = 2 * n1 - 2;

        return TestResults<T>(t, df, statsTest(t, TDist<T>(df), type));
    }

    /*
     * Paired two-groups t-test
     */

    template <typename T, typename Iter> static TestResults<T> tt_paired_two_groups
            (const Iter &x, const Iter &y, TestType type = TwoTailed, T h0 = T())
    {
        std::vector<T> diffs;
        diffs.resize(x.size());
        
        for (auto i = 0; i < x.size(); i++)
        {
            diffs[i] = y[i] - x[i];
        }
        
        const auto s  = sd(diffs);
        const auto df = x.size() - 1;
        const auto t  = s / sqrt(x.size());

        return TestResults<T>(t, df, statsTest(t, TDist<T>(df), type));
    }

    /*
     * Welch's t-test. An adaptation of Student's t-test and is more robust when the two
     * samples have unequal variances and unequal sample sizes (heteroscedastic).
     */

    template <typename T, typename Iter> static TestResults<T> tt_welch
        (const Iter &x, const Iter &y, TestType type = TwoTailed, T h0 = T())
    {
        const auto u1 = mean(x);
        const auto u2 = mean(y);
        const auto s1 = sd(x);
        const auto s2 = sd(y);
        const auto n1 = x.size();
        const auto n2 = y.size();
        const auto s  = std::sqrt((s1 * s1 / n1) + (s2 * s2 / n2));
        const auto t  = (u1 - u2 - h0) / s;

        const DF df = std::pow(s1 * s1 / n1 + s2 * s2 / n2, 2) /
                                (std::pow(s1 * s1 / n1, 2) / (n1 - 1) +
                                 std::pow(s2 * s2 / n2, 2) / (n2 - 1));

        return TestResults<T>(t, df, statsTest(t, TDist<T>(df), type));
    }

    template <typename T> struct TTest
    {
        template <typename Iter> static TestResults<T> test
            (const Iter &x, const Iter &y, TestType type = TwoTailed, T h0 = T())
        {
            return tt_welch(x, y, type, h0);
        }
    };

    static TTest<double> t;
}

#endif