#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

#include <numeric>
#include "data/tokens.hpp"
#include "MetaQuin/m_blat.hpp"
#include "stats/analyzer.hpp"
#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAssembly
    {
        enum Software
        {
            MetaQuast,
        };

        struct Stats : public LinearStats
        {
            // Statistics for de-novo assembly
            DAsssembly::Stats<DAsssembly::Contig> dnovo;
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            // Eg: alignments_contigs.tsv (MetaQuast)
            FileName contigs;
            
            // Eg: genome_info.txt (MetaQuast)
            FileName genome;

            Software soft;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif