#ifndef GI_F_VIEWER_HPP
#define GI_F_VIEWER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FViewer
    {
        struct Options : public ViewerOptions
        {
            // Empty Implementation
        };

        static void report(const std::string &, const ViewerOptions &options = ViewerOptions());
    };
}

#endif