#ifndef GI_F_VIEWER_HPP
#define GI_F_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FViewer
    {
        typedef ViewerOptions Options;
        
        // Generate a IGV session for fusion analysis
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("Fusion " + o.path + " " + file);
        }
    };
}

#endif