#ifndef SS_MATRIX_HPP
#define SS_MATRIX_HPP

#include <Eigen/Dense>
#include <ss/data/data.hpp>

/*
 * This file defines convenient functions for manipulating an Eigen matrix.
 */

namespace SS
{
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;
    
    namespace MatrixUtils
    {
        // Apply a functor to every element in the matrix
        template <typename F> void each(const Matrix &m, F f)
        {
            for (auto i = 0; i < m.rows(); i++)
            {
                for (auto j = 0; j < m.cols(); j++)
                {
                    f(i, j, m(i, j));
                }
            }
        }

        template <typename Iter> static Vector vector(const Iter &x)
        {
            Vector X(x.size());
            
            auto i = 0;
            
            for (const auto &t : x)
            {
                X(i++) = t;
            }

            return X;
        }
        
        template <typename T> static Matrix matrix(const T *t, Counts n)
        {
            Matrix X(n, 1);

            for (auto i = 0; i < n; i++, t++)
            {
                X(i, 0) = *t;
            }
            
            return X;
        }
        
        template <typename Iter> static Matrix matrix(const Iter &x)
        {
            Matrix X(x.size(), 1);
            
            auto i = 0;
            
            for (const auto &t : x)
            {
                X(i++, 0) = t;
            }
            
            return X;
        }

        struct MatrixCount
        {
            // Row sums
            Eigen::VectorXd rsums;
            
            // Column sums
            Eigen::VectorXd csums;

            // Overall sums
            Real sums = 0;
        };

        inline MatrixCount count(const Matrix &m)
        {
            MatrixCount c;
            
            c.rsums = Eigen::VectorXd(m.rows());
            c.csums = Eigen::VectorXd(m.cols());

            c.rsums.setConstant(0);
            c.csums.setConstant(0);
            
            MatrixUtils::each(m, [&](Index i, Index j, Real x)
            {
                c.sums += x;
                c.rsums(i) += x;
                c.csums(j) += x;
            });
            
            return c;
        }
    }
}

#endif