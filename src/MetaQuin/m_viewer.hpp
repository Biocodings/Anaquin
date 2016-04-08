#ifndef M_VIEWER_HPP
#define M_VIEWER_HPP

#include "data/script.hpp"

namespace Anaquin
{
    typedef ViewerOptions Options;
    
    static void generate(const FileName &file, const Options &o = Options())
    {
        Script::viewer("MetaQuin " + o.path + " " + file);
    }
}

#endif