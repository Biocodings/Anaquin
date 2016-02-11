#ifndef F_VIEWER_HPP
#define F_VIEWER_HPP

#include "data/script.hpp"

namespace Anaquin
{
    struct FViewer
    {
        typedef ViewerOptions Options;
        
        static void generate(const FileName &file, const Options &o = Options())
        {
            Script::viewer("Fusion " + o.path + " " + file);
        }
    };
}

#endif