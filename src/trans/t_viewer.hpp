#ifndef GI_T_VIEWER_HPP
#define GI_T_VIEWER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TViewer
    {
        typedef ViewerOptions Options ;

        static void generate(const std::string &, const Options &options = Options());
    };
}

#endif