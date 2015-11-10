#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

#include <numeric>
#include "data/tokens.hpp"
#include "meta/m_blat.hpp"
#include "stats/analyzer.hpp"
#include "meta/m_assembler.hpp"

namespace Anaquin
{
    enum MetaAssembler
    {
        Velvet,
        RayMeta,
    };
    
    struct MAssembly
    {
        struct Stats : public DAsssembly::Stats<Contig>
        {
            // Statistics for the alignment
            MBlat::Stats blat;
            
            // Statistics for abundance
            LinearStats lm;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            /*
             * Metagenomics assembler report results differently. For example, RayMeta generates
             * "Contigs.tsv" which specifies the coverage of the contigs.
             */
            
            FileName contigs;
            
            // Alignment file by blat
            FileName psl;

            // The type of the assembler used
            MetaAssembler tool = MetaAssembler::Velvet;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif