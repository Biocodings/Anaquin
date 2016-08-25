#ifndef SS_ANOVA_HPP
#define SS_ANOVA_HPP

#include <ss/data/results.hpp>
#include <ss/data/variable.hpp>
#include <ss/internal/model.hpp>
#include <ss/regression/linear.hpp>

namespace SS
{
    /*
     *                      Analysis of Variance (ANOVA)
     *
     * Analysis of variance (ANOVA) is a collection of statistical models used to analyze the
     * differences among group means and their associated procedures. In the ANOVA setting, the
     * observed variance in a particular variable is partitioned into components attributable to
     * different sources of variation.
     */

    template <typename... Args> StatsResults anova(const Variable &x, Args... args)
    {
        StatsResults r("anova");

        Internal::vArgs([&](std::vector<const Variable *> &p)
        {
            // Dependent variable
            const auto Y = MatrixUtils::matrix(p[0]->data(), p[0]->size());
  
            p.erase(p.begin()+0);
            
            // Independnet variables
            const auto X = Internal::modelMatrix(p);

            // Perform a regression on the dummy variables
            const auto lm = linearModel(Y, X);

            using namespace Keys;

            r[FStats]  = lm.f;
            r[PValue]  = lm.p;
            r[ModelDF] = lm.model.df;
            r[ModelSS] = lm.model.ss;
            r[ModelMS] = lm.model.ms;
            r[ErrorDF] = lm.error.df;
            r[ErrorSS] = lm.error.ss;
            r[ErrorMS] = lm.error.ms;
            r[TotalDF] = lm.total.df;
            r[TotalSS] = lm.total.ss;
            r[TotalMS] = lm.total.ms;

        }, x, args...);

        return r;
    }
}

#endif