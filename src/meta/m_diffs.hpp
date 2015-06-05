#ifndef GI_M_DIFFS_HPP
#define GI_M_DIFFS_HPP

#include "analyzer.hpp"

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
            // Empty Implementation
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif