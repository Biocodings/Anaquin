#ifndef GI_T_VIEWER_HPP
#define GI_T_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TViewer
    {
        typedef ViewerOptions Options;

        // Generate a IGV session for transcriptome analysis
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("Trans " + o.path + " " + file);
        }
    };
}

#endif