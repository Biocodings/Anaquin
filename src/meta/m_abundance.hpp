#ifndef GI_M_ABUNDANCE_HPP
#define GI_M_ABUNDANCE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAbundance
    {
        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif