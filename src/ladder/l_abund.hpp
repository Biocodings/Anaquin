#ifndef GI_L_ABUND_HPP
#define GI_L_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LAbund
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : ModelStats
        {
            // Histogram expected
            std::map<TypeID, Coverage> expect;

            // Histogram before normalization and adjustment
            std::map<TypeID, Counts> abund;

            // Histogram after normalization but before adjustment
            std::map<TypeID, Coverage> actual;

            // Histogram after adjustment
            std::map<TypeID, Coverage> adjusted;
            
            // Adjusted abundance for each sequin
            std::map<SequinID, Coverage> s_adjusted;

            SequinCounter c;
            
            // Expected size of the library
            Counts expTotal = 0;

            // Measured size of the library
            Counts actTotal = 0;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif