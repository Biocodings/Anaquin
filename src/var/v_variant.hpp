#ifndef GI_V_VARIANT_HPP
#define GI_V_VARIANT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VVariant
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            // Overall performance
            Confusion m;

            // The proportion of variations with alignment coverage
            double covered;

            // Measure of variant detection independent to sequencing depth or coverage
            double efficiency;

            SequinHist c = Analyzer::histogram(Standard::instance().v_seqs_A);
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif