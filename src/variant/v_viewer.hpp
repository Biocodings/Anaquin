#ifndef GI_V_VIEWER_HPP
#define GI_V_VIEWER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VViewer
    {
        typedef SingleMixtureOptions Options;
        
        struct Stats
        {
            // Empty Implementation
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif