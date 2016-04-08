#ifndef V_VIEWER_HPP
#define V_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VViewer
    {
        typedef ViewerOptions Options;
        
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("VarQuin " + o.path + " " + file);
        }
    };
}

#endif