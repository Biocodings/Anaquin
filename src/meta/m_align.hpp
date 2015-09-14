#ifndef GI_M_ALIGN_HPP
#define GI_M_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAlign
    {
        static Stats report(const std::string &file, const Options &options = Options());
    };
}

#endif