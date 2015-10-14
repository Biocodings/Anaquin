#ifndef F_ALIGN_HPP
#define F_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class FAlign : public Analyzer
    {
        public:

            struct Stats : public AlignmentStats
            {
            };
        
            typedef AnalyzerOptions Options;

            static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif