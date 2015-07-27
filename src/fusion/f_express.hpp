#ifndef GI_F_EXPRESS_HPP
#define GI_F_EXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FExpress
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        static ModelStats analyze(const std::string &, const Options &options = Options());
    };
}

#endif