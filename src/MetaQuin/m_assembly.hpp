#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

#include <numeric>
#include "data/tokens.hpp"
#include "MetaQuin/m_blat.hpp"
#include "stats/analyzer.hpp"
#include "MetaQuin/m_assembler.hpp"

namespace Anaquin
{
    struct MAssembly
    {
        enum Software
        {
            MetaQuast,
        };

        struct Stats : public DAsssembly::Stats<Contig>
        {
            // Statistics for the alignment
            MBlat::Stats blat;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            Software soft;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif