#ifndef GI_C_JOIN_HPP
#define GI_C_JOIN_HPP

#include "stats/analyzer.hpp"

namespace Spike
{
    struct CJoin
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
            std::map<SequinID, Counts> raw;

            // Histogram after normalization but before correction
            std::map<SequinID, Coverage> normal;
            
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