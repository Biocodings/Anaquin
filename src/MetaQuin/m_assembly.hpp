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
            Velvet,
            RayMeta,
        };

        struct Stats : public DAsssembly::Stats<Contig>
        {
            // Statistics for the alignment
            MBlat::Stats blat;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            /*
             * The assemblers report results differently. For example, RayMeta generates
             * "Contigs.tsv" that specifies the coverage of the contigs. This is defined
             * for RayMeta.
             */
            
            FileName contigs;
            
            // Alignment file by blat
            FileName psl;

            Software soft;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif