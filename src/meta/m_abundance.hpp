#ifndef GI_M_ABUNDANCE_HPP
#define GI_M_ABUNDANCE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAbundance
    {
        struct Stats
        {
            
        };
        
        struct Options : public AnalyzerOptions
        {

        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif