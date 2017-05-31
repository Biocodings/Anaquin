#ifndef V_CONJOIN_HPP
#define V_CONJOIN_HPP

#include <map>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VConjoint
    {
        typedef AnalyzerOptions Options;
        
        struct Stats
        {
            struct Data
            {
                Counts  measured = 0;
                Concent expected;
            };

            std::map<SequinID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
