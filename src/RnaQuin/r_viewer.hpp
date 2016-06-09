#ifndef R_VIEWER_HPP
#define R_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TViewer
    {
        typedef ViewerOptions Options;

        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("TransQuin " + o.path + " " + file);
        }
    };
}

#endif