#ifndef GI_V_DISCOVER_HPP
#define GI_V_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        typedef FuzzyOptions Options;

        struct Stats : public MappingStats
        {
            long detected = 0;
            
            // Overall performance
            Confusion m;

            SequinHist h = Standard::instance().r_var.hist();
        };

        static Stats report(const std::string &, const Options &options = Options());
    };
}

#endif