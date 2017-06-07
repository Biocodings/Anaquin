#ifndef V_CONJOINT_HPP
#define V_CONJOINT_HPP

#include <map>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VConjoint
    {
        typedef AnalyzerOptions Options;

        struct Stats
        {
            std::map<SequinID, Counts> data;
        };

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
