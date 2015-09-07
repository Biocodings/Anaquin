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
            long detected;
            
            // Overall performance
            Confusion m;

            // The proportion of variations with alignment coverage
            double covered;

            // Measure of variant detection independent to sequencing depth or coverage
            double efficiency;

            SequinHist h = Standard::instance().r_var.hist();
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif