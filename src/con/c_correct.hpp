#ifndef GI_C_CORRECT_HPP
#define GI_C_CORRECT_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct CCorrect
    {
        struct Options : public SingleMixtureOptions
        {
            // Empty Implementation
        };

        struct Stats : public ModelStats
        {
            // Histogram expected
            std::map<SequinID, Coverage> expect;

            // Histogram before normalization and correction
            std::map<SequinID, Counts> abund;

            // Histogram after normalization but before correction
            std::map<SequinID, Coverage> actual;

            // Histogram after correction
            std::map<SequinID, Coverage> correct;
            
            // Expected size of the library
            Counts expTotal = 0;

            // Measured size of the library
            Counts actTotal = 0;
        };

        static Stats analyze(const std::string &file, const Options &options = Options());
    };
}

#endif