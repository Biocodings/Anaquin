#ifndef GI_L_CORRECT_HPP
#define GI_L_CORRECT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LCorrect
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats
        {
            // Histogram expected
            std::map<TypeID, Coverage> expect;

            // Histogram before normalization and correction
            std::map<TypeID, Counts> abund;

            // Histogram after normalization but before correction
            std::map<TypeID, Coverage> actual;

            // Histogram after correction
            std::map<TypeID, Coverage> correct;
            
            // Corrected abundance for each sequin
            std::map<SequinID, Coverage> s_correct;
            
            // Expected size of the library
            Counts expTotal = 0;

            // Measured size of the library
            Counts actTotal = 0;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif