#ifndef T_KEXPRESS_HPP
#define T_KEXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TKExpress : public Analyzer
    {
        enum class Software
        {
            Kallisto,
        };
        
        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            // Gene or isoform?
            Metrics metrs;
            
            // Only Kallisto is supported
            Software soft;
        };

        struct Stats : public MappingStats, public SequinStats, public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o);
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
