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

        struct Stats : public SequinStats
        {
            Limit limit;
        };

        static Scripts generateCSV(const Stats &, const Options &);

        static Scripts generateSummary(const FileName &,
                                       const Stats &,
                                       const Options &);

        static Scripts generateRLinear(const FileName &,
                                       const Stats &,
                                       const Options &);

        static void writeSummary(const FileName &,
                                 const FileName &,
                                 const Stats &,
                                 const Options &);
        
        static void writeCSV(const FileName &, const Stats &, const Options &);
        static void writeRLinear(const FileName &, const FileName &, const Stats &, const Options &);

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
