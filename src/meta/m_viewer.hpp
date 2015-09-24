#ifndef M_VIEWER_HPP
#define M_VIEWER_HPP

#include "data/script.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    typedef ViewerOptions Options;
    
    // Generate a IGV session for metagenomcis analysis
    static void generate(const FileName &file, const Options &o = Options())
    {
        Script::viewer("Meta " + o.path + " " + file);
    }
}

#endif