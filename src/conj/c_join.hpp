#ifndef GI_C_JOIN_HPP
#define GI_C_JOIN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct CJoin
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            std::map<SequinID, Counts> hist;
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif