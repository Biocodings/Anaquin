#ifndef V_STRUCTURE_HPP
#define V_STRUCTURE_HPP

#include "data/vData.hpp"
#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VStructure
    {
        typedef AnalyzerOptions Options;
        
        struct EStats
        {
            unsigned found = 0;
        };
        
        struct SStats
        {
        };

        static EStats analyzeE(const FileName &, const Options &o);
        static SStats analyzeS(const FileName &, const Options &o);

        // Report for both endogenous and sequins
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
