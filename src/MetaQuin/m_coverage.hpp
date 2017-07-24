#ifndef M_COVERAGE_HPP
#define M_COVERAGE_HPP

#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MCoverage
    {
        struct Stats : public SequinStats, public AlignmentStats, public HistStats
        {
            // Empty Implementation
        };
        
        enum class Format
        {
            BAM,
            RayMeta,
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            Format format;
            Mixture mix = Mixture::Mix_1;
        };

        static Stats analyze(const std::vector<FileName> &, const Options &);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
