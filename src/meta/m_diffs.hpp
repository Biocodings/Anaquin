#ifndef GI_M_DIFFS_HPP
#define GI_M_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct MDiffs
    {
        struct Stats : public CorrelationStats
        {
            // Empty Implementation
        };

        struct Options : public DoubleMixtureOptions
        {
            // An optional PSL file for the first sample
            std::string psl_1;

            // An optional PSL file for the second sample
            std::string psl_2;
        };

        static Stats analyze(const std::string &, const std::string &, const Options &options = Options());
    };
}

#endif