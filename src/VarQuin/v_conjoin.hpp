#ifndef V_CONJOIN_HPP
#define V_CONJOIN_HPP

#include "data/vData.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VConjoint
    {
        typedef AnalyzerOptions Options;
        
        struct Stats
        {

        };

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
