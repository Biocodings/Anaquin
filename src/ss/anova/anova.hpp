#ifndef SS_ANOVA_HPP
#define SS_ANOVA_HPP

#include <ss/r.hpp>
#include <ss/anova/coding.hpp>
#include <ss/regression/lm.hpp>

namespace SS
{
    template<typename T> struct AnovaResults
    {
        // Total variation
        Variation<T> total;
        
        // Variation by the residuals
        Variation<T> error;
        
        // Variation by the model
        Variation<T> model;
        
        P p;
        T f;
    };
    
    /*
     * Fit an analysis of variance model by a call to lm for each stratum.
     *
     * This provides a wrapper to lm for fitting linear models to balanced or unbalanced
     * experimental designs.
     */
    
    template<typename T, typename Coding> AnovaResults<T> aov(const Formula &f, const Internal::Data<T> &d)
    {
        const auto X = Coding().encode(d.cs[0]._data, d.cs[1]._data);
        const auto Y = convert(std::vector<Internal::C<T>> { d.cs[0] }, Dependent);
        
        // Wrapping the call to its equivalent linear model
        const auto rm = linearModel(f, *Y, X, mean(d.cs[0]._data));

        AnovaResults<T> r;
        
        r.model.df = rm.model.df;
        r.model.ss = rm.model.ss;
        r.model.ms = rm.model.ms;
        
        r.error.df = rm.error.df;
        r.error.ss = rm.error.ss;
        r.error.ms = rm.error.ms;
        
        r.total.df = rm.total.df;
        r.total.ss = rm.total.ss;
        r.total.ms = rm.total.ms;
        
        r.f = rm.f;
        r.p = rm.p;
        
        return r;
    }

    template <typename Coding = DummyCoding> AnovaResults<Real> aov(const Formula &f, const Internal::Data<Real> &d)
    {
        return aov<Real, Coding>(f, d);
    }
}

#endif