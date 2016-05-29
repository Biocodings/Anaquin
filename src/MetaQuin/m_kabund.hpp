#ifndef M_KABUND_HPP
#define M_KABUND_HPP

#include "MetaQuin/m_assembly.hpp"

namespace Anaquin
{
    struct MKAbund
    {   
        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            MBlat::Stats blat;
            
            // Statistics for the assembly
            DAsssembly::Stats<DAsssembly::Contig> assembly;
        };

        struct Options : public AnalyzerOptions
        {
            // Required by the GCC compiler ...
            Options() {}

            MAligner aligner;
            MAssembler assembler;
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &o = Options());
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif