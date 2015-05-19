#ifndef SS_MATRIX_HPP
#define SS_MATRIX_HPP

#include <vector>
#include <Eigen/Dense>
#include <ss/types.hpp>

namespace SS
{
    typedef Eigen::MatrixXd Matrix;

    template <typename F> void for_each(const Matrix &m, F f)
    {
        for (auto i = 0; i < m.rows(); i++)
        {
            for (auto j = 0; j < m.cols(); j++)
            {
                f(i, j, m(i, j));
            }
        }
    }

    template<typename T> struct MatrixCount
    {
        // Row sums
        std::vector<T> rs;
        
        // Column sums
        std::vector<T> cs;

        // Overall sums
        T sums = T();
    };

    inline void count(const Matrix &m, MatrixCount<Real> &c)
    {
        c.rs.clear();
        c.cs.clear();
        c.rs.resize(m.rows());
        c.cs.resize(m.cols());

        for_each(m, [&](std::size_t i, std::size_t j, Real x)
        {
            c.sums  += x;
            c.rs[i] += x;
            c.cs[j] += x;
        });
    }
}

#endif