#ifndef GI_V_VIEWER_HPP
#define GI_V_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VViewer
    {
        typedef ViewerOptions Options;
        
        // Generate a IGV session for variant analysis
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("Variant " + o.path + " " + file);
        }
    };
}

#endif