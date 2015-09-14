#ifndef GI_M_ABUNDANCE_HPP
#define GI_M_ABUNDANCE_HPP

#include "meta/m_assembly.hpp"

namespace Anaquin
{
    struct MAbundance
    {
        typedef MAssembly::Stats   Stats;
        typedef MAssembly::Options Options;

        static Stats report(const std::string &, const Options &o = Options());
    };
}

#endif