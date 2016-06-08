#ifndef V_KALLELE_HPP
#define V_KALLELE_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_freq.hpp"

namespace Anaquin
{
    struct VKAllele
    {
        typedef IndexOptions Options;
        typedef VFreq::Stats Stats;

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif