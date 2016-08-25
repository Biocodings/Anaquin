#ifndef SS_SEGMENTED_HPP
#define SS_SEGMENTED_HPP

#include <map>
#include <limits>
#include <ss/data/errors.hpp>
#include <ss/internal/rank.hpp>
#include <ss/regression/linear.hpp>

namespace SS
{
    template <typename T> struct SegmentedResults
    {
        // Sums of residuals for both regression
        std::map<T, T> sums;

        // Correlation for both regressions
        std::map<T, T> lr, rr;

        // R2 for both regressions
        std::map<T, T> lR2, rR2;

        // Slope for both regressions
        std::map<T, T> lSl, rSl;
        
        // Intercept for both regressions
        std::map<T, T> lInt, rInt;

        // The break point
        T b = NAN;
        
        // Index of the break point
        unsigned index;

        // Left correlation at the break-point
        inline T blr() const { return lr.at(b); }
        
        // Right correction at the break-point
        inline T brr() const { return rr.at(b); }
        
        // Left R2 at the break-point
        inline T blR2() const { return lR2.at(b); }
            
        // Right R2 at the break-point
        inline T brR2() const { return rR2.at(b); }

        // Left slope at the break-point
        inline T blSl() const { return lSl.at(b); }
        
        // Right slope at the break-point
        inline T brSl() const { return rSl.at(b); }
        
        // Left intercept at the break-point
        inline T blInt() const { return lInt.at(b); }
        
        // Right intercept at the break-point
        inline T brInt() const { return rInt.at(b); }
    };

    template <typename T> SegmentedResults<T> segmentPearson(const std::vector<T> &x, const std::vector<T> &y)
    {
        SS_ASSERT(!x.empty() && !y.empty(), "Failed to fit a segmented regression on empty inputs");
        
        SegmentedResults<T> r;
        
        auto xs = x;
        auto ys = y;
        
        Internal::permSort(xs, ys);
        
        auto max = std::numeric_limits<double>::min();
        
        for (auto i = 3; i < x.size() - 3; i++)
        {
            const std::vector<T> x1(xs.begin(), xs.begin()+i);
            const std::vector<T> y1(ys.begin(), ys.begin()+i);
            const std::vector<T> x2(xs.begin()+i, xs.end());
            const std::vector<T> y2(ys.begin()+i, ys.end());
            
            assert(x1.size() + x2.size() == x.size());
            assert(y1.size() + y2.size() == y.size());
            
            const auto m1 = linearModel(y1, x1);
            const auto m2 = linearModel(y2, x2);
            
            // Breakpoint for iteration
            const auto b = xs[i-1];
            
            r.lr[b] = corrPearson(x1, y1);
            r.rr[b] = corrPearson(x2, y2);
            
            // Add up the residuals
            r.sums[b] = m1.error.ss + m2.error.ss;
            
            r.lR2[b] = m1.r2;
            r.rR2[b] = m2.r2;
            
            r.lSl[b] = m1.coeffs[1].est;
            r.rSl[b] = m2.coeffs[1].est;
            
            r.lInt[b] = m1.coeffs[0].est;
            r.rInt[b] = m2.coeffs[0].est;
            
            if (r.rr[b] > max)
            {
                max = r.sums[r.b = b];
            }
        }
        
        return r;
    }

    template <typename T> SegmentedResults<T> segmentSSE(const std::vector<T> &x, const std::vector<T> &y)
    {
        SS_ASSERT(!x.empty() && !y.empty(), "Failed to fit a segmented regression on empty inputs");

        SegmentedResults<T> r;
        
        auto xs = x;
        auto ys = y;
        
        Internal::permSort(xs, ys);
        
        /*
         * Estimate all possible regression models
         */
        
        auto min = std::numeric_limits<double>::max();
        auto j = 0;

        for (auto i = 3; i < x.size() - 3; i++)
        {
            const std::vector<T> x1(xs.begin(), xs.begin()+i);
            const std::vector<T> y1(ys.begin(), ys.begin()+i);
            const std::vector<T> x2(xs.begin()+i, xs.end());
            const std::vector<T> y2(ys.begin()+i, ys.end());

            assert(x1.size() + x2.size() == x.size());
            assert(y1.size() + y2.size() == y.size());

            // Breakpoint for iteration
            const auto b = xs[i-1];
            
            // Only estimate a model once for any x-value
            if (r.sums.count(b))
            {
                continue;
            }

            const auto m1 = linearModel(y1, x1);
            const auto m2 = linearModel(y2, x2);

            r.lr[b] = corrPearson(x1, y1);
            r.rr[b] = corrPearson(x2, y2);
            
            // Add up the residuals
            r.sums[b] = m1.error.ss + m2.error.ss;

            r.lR2[b] = m1.r2;
            r.rR2[b] = m2.r2;
            
            r.lSl[b] = m1.coeffs[1].est;
            r.rSl[b] = m2.coeffs[1].est;

            r.lInt[b] = m1.coeffs[0].est;
            r.rInt[b] = m2.coeffs[0].est;
            
            if (r.sums[b] < min)
            {
                min = r.sums[r.b = b];
                r.index = j;
            }
            
            j++;
        }

        return r;
    }
}

#endif