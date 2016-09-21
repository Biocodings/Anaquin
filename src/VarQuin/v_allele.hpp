#ifndef V_ALLELE_HPP
#define V_ALLELE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAllele
    {
        enum class Format
        {
            Salmon
        };
        
        struct Options : public AnalyzerOptions
        {
            Format format;
        };

        struct Stats : public LinearStats
        {
            // Empty Implementation
        };

        static void generateCSV(const FileName &, const Stats &, const Options &);
        static void generateSummary(const FileName &, const Stats &, const Options &);
        static void generateR(const FileName &, const FileName &, const Stats &, const Options &);

        static Stats analyze(const FileName &, const Options &o);
        static void report  (const FileName &, const Options &o = Options());
    };
}

#endif