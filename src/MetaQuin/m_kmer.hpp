#ifndef M_KMER_HPP
#define M_KMER_HPP

#include "MetaQuin/m_assembly.hpp"

namespace Anaquin
{
    struct MKMer
    {
        enum class Software
        {
            Kallisto
        };
        
        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            // Empty Implementation
        };

        struct Options : public AnalyzerOptions
        {
            // Required by the GCC compiler ...
            Options() {}

            Software soft;
        };
        
        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif