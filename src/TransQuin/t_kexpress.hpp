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
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            // Only Kallisto is supported
            Software soft;
        };

        struct Stats : public MappingStats, public SequinStats, public LinearStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o);
        
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }

            return stats;
        }
        
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
