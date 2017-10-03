#ifndef V_PROCESS_HPP
#define V_PROCESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VProcess
    {
        typedef AnalyzerOptions Options;

        enum class Status
        {
            ReverseReverse,
            ReverseNotMapped,
            ForwardForward,
            ForwardReverse,
            ForwardNotMapped,
            NotMappedNotMapped,
            RevHang,
            ForHang
        };

        struct Stats : public MappingStats
        {
            std::map<VProcess::Status, Counts> counts;
        };

        static Stats analyze(const FileName &, const Options &);
        static void report(const FileName &, const Options &);
    };
}

#endif
