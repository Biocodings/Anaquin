#ifndef GI_T_VIEWER_HPP
#define GI_T_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TViewer
    {
        typedef ViewerOptions Options;

        static void generate(const std::string &, const Options &o = Options())
        {
            Script::viewer("Trans " + o.path + " " + o.file);
        }
    };
}

#endif