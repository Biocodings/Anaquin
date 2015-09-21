#ifndef GI_M_ASSEMBLY_HPP
#define GI_M_ASSEMBLY_HPP

#include <numeric>
#include "data/tokens.hpp"
#include "meta/m_blat.hpp"
#include "meta/m_velvet.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MAssembly
    {
        enum Assembler
        {
            Velvet,
        };

        struct Stats : public DAsssembly::Stats<Contig>
        {
            // Statistics for the alignment
            MBlat::Stats blat;
            
            // Statistics for abundance
            LinearStats lm;
        };

        struct Options : public SingleMixtureOptions
        {
            Options() {}

            // Alignment file by blat
            FileName psl;

            // The type of the assembler used
            Assembler tool = Assembler::Velvet;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static Stats report (const FileName &, const Options &o = Options());
    };
}

#endif