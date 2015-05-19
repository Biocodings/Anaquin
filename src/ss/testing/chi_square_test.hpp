#ifndef SS_CHI_SQUARE_TEST
#define SS_CHI_SQUARE_TEST

#include <ss/matrix.hpp>
#include <ss/testing/test.hpp>
#include <ss/dists/chi_square.hpp>

namespace SS
{
    template<typename T> TestResults<T> independence(const Matrix &m)
    {
        MatrixCount<T> c;
        
        // We need the marginal distribution for the test statistic
        count(m, c);
        
        TestStats x = TestStats();
        
        for_each(m, [&](std::size_t i, std::size_t j, const T &o)
        {
            /*
             * The expected frequency under null hypothesis is the product of expected marginal distributions.
             */
            
            const auto e = (c.rs[i] * c.cs[j] / c.sums);
            
            x += std::pow(e - o, 2) / e;
        });

        return stats_test<TestResults<T>, Internal::ChiSquare<T>>(x, ChiSquare(c.rs.size() * c.cs.size() - 1), TwoTailed);
    }
}

#endif