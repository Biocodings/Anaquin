#ifndef L_IGV_HPP
#define L_IGV_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LViewer
    {
        typedef ViewerOptions Options;
        
        // Generate a IGV session for ladder analysis
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("Ladder " + o.path + " " + file);
        }
    };
}

#endif