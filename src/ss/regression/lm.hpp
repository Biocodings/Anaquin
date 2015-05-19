#ifndef SS_LM_HPP
#define SS_LM_HPP

#include <ss/data.hpp>
#include <ss/matrix.hpp>
#include <ss/dists/f_dist.hpp>
#include <ss/testing/t_test.hpp>

namespace SS
{
    typedef std::string Formula;
    
    template <typename T> struct Variation
    {
        // The sum-of-squares of the variation
        T ss;
        
        // Mean square of the variation
        T ms;
        
        DF df;
    };
    
    template <typename T> struct Coefficient
    {
        T v;
        
        // Test statistic whether the coefficient is significant from zero
        T t;
        
        // Standard error of the coefficient
        T e;
        
        // P-value for the test statistic
        P p;
        
        DF dt;
    };
    
    template <typename T> struct RegressionResults
    {
        std::vector<Coefficient<T>> coeffs;
        
        Variation<T> error;
        Variation<T> model;
        Variation<T> total;
        
        Percentage r2;
        Percentage ar2;
        
        // The F-statistic for the regression
        T f;
        
        // The p-value for the F-statistic
        P p;
    };
    
    enum InputType
    {
        Independent,
        Dependent
    };
    
    /*
     * Convert a given input to its matrix representation. The exact representation depends
     * on the conversion type. For example
     */
    
    template <typename T> std::shared_ptr<Matrix> convert(const std::vector<Internal::C<T>> &cs, InputType type)
    {
        assert(!cs.empty());
        
        switch (type)
        {
            case Independent:
            {
                std::shared_ptr<Matrix> X(new Matrix(cs[0].size(), cs.size() + 1));
                
                // Constant coefficient
                X->col(0).setOnes();
                
                for (auto i = 0; i < cs.size(); i++)
                {
                    for (auto j = 0; j < cs[i].size(); j++)
                    {
                        (*X)(j, i+1) = cs[i]._data[j];
                    }
                }
                
                return X;
            }
                
            case Dependent:
            {
                std::shared_ptr<Matrix> Y(new Matrix(cs[0].size(), 1));
                
                for (auto i = 0; i < cs[0].size(); i++)
                {
                    (*Y)(i, 0) = cs[0]._data[i];
                }
                
                return Y;
            }
        }
    }
    
    template <typename T> RegressionResults<T> linearModel(const Formula &f, const Matrix &Y, const Matrix &X, T mean)
    {
        assert(Y.rows() == X.rows());
        using namespace Eigen;
        
#define DEFINE_EY() Matrix EY(Y.rows(), 1); EY.setConstant(mean);
        DEFINE_EY();
        
        const Matrix X_ = X.transpose();
        const Matrix H = (X_ * X);
        
        // Covariance matrix (multiple by sigma^2 gives the covariance)
        const Matrix C = H.inverse();
        
        // Fitted coefficients
        const Matrix B = (C * X_) * Y;
        
        // Residuals from the the model (Nx1 column vector)
        const Matrix E = Y - X * B;
        
        RegressionResults<T> r;
        
        // The size of the sample
        const auto n = Y.rows();
        
        // Number of coefficients (including the constant coefficient)
        const auto p = X.cols();
        
        /*
         * Sum of squares of error (SSE)
         *
         *    SSE = E'E
         */
        
        r.error.df = n - p;
        r.error.ss = (E.transpose() * E)(0, 0);
        r.error.ms = r.error.df ? r.error.ss / r.error.df : NAN;
        
        /*
         * Sum of squares of model (SSM)
         *
         *    SSM = (XB-EY)'(XB-EY)
         */
        
        const Matrix M = (X * B) - EY;
        
        r.model.df = p - 1;
        r.model.ss = (M.transpose() * M)(0, 0);
        r.model.ms = r.model.ss / r.model.df;
        
        assert(r.model.ss >= 0);
        assert(r.model.ms >= 0);
        
        r.total.df = r.model.df + r.error.df;
        r.total.ss = r.model.ss + r.error.ss;
        r.total.ms = r.total.ss / r.total.df;
        
        /*
         * Compute the F-statistic for the regression. It tests for a significant linear regression
         * relationship between the independent variable and the predictor variables.
         */
        
        r.f = r.model.ms / r.error.ms;
        r.p = P(r.error.df ? 1 - FDist<T>(r.model.df, r.error.df).cdf(r.f) : NAN);
        r.r2  = r.model.ss / r.total.ss;
        r.ar2 = 1 - (r.error.ss / (n - p)) / (r.total.ss / (n - 1));
        
        /*
         * Variance of the coefficients. It's a product of the gaussian noise and
         * the distance of the indepndent variables to their average.
         */
        
        const auto V = r.error.ms * C;
        
        /*
         * Compute statistics for each coefficient. Each coefficient in the model is a Gaussian
         * random variable. The magnitude is the estimate of the mean of the distribution, and
         * the standard error is the square root of the variance of that distribution.
         */
        
        for (auto i = 0; i < p; i++)
        {
            Coefficient<T> c;
            
            /*
             * Under the assumption of normality of the error terms, the estimator of the slope
             * coefficient will itself be normally distributed with mean B(i,0) and variance V(i,i).
             */
            
            c.v = B(i, 0);
            c.e = sqrt(V(i, i));
            
            // Each coefficient inherits the same degree of freedom from MSE
            c.dt = r.error.df;
            
            /*
             * Under the null hypothesis, each coefficent has a value of zero. The coefficents
             * are normally distributed while the population coefficient is unknown.
             */
            
            c.t = c.v / c.e;
            c.p = c.dt ? statsTest(c.t, TDist<T>(c.dt), TwoTailed).p : NAN;
            
            r.coeffs.push_back(c);
        }
        
        return r;
    }
    
    template <typename T> RegressionResults<T> lm(const Formula &f, const Internal::Data<T> &d)
    {
        const auto &y = d.cs[0]._data;
        const auto &x = d.cs[1]._data;
        
        if (y.size() != x.size())
        {
            throw std::runtime_error("Mismatch size in samples");
        }
        
        const auto X = convert(std::vector<Internal::C<T>> { d.cs[1] }, Independent);
        
        #define DEFINE_Y()  Matrix Y(y.size(), 1); for (auto i = 0; i < y.size(); i++) { Y(i, 0) = y[i]; }
        DEFINE_Y();
        
        return linearModel(f, Y, *X, mean(y));
    }
}

#endif