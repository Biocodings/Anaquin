#ifndef M_ABUND_HPP
#define M_ABUND_HPP

#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAbund
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
            Format format;
        };
        
        static Scripts generateRLinear(const FileName &, const Stats &, const Options &);
        
        static Stats analyze(const std::vector<FileName> &, const Options &);
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
