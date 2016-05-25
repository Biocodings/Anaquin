#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAssembly
    {
        enum Assembler
        {
            Velvet,
            RayMeta,
        };

        enum Alignment
        {
            MetaQuast,
            Blat,
        };

        struct Stats : public LinearStats, public DAsssembly::DenoAssemblyImpl
        {
            inline SequinID findSeq(const ContigID &cID) const
            {
                auto x = cID;
                
                if (soft == RayMeta && align == Blat)
                {
                    /*
                     * BLAT doesn't do spacing. For example, "contig-0 52976 nucleotides"
                     * would be given as "contig-0".
                     */
                    
                    x = Tokens::first(cID, " ");
                }
                
                return c2s.count(x) ? c2s.at(x) : "";
            }
            
            /*
             * Usually, they're not needed but we need them here to check for special cases.
             * Commented in findSeq().
             */

            Assembler soft;
            Alignment align;

            std::map<ContigID, SequinID> c2s;
            std::map<SequinID, std::vector<ContigID>> s2c;
            
            // Statistics for de-novo assembly
            DAsssembly::Stats<DAsssembly::Contig> dnovo;
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Assembler soft;
            Alignment align;
        };

        static Stats analyze(const std::vector<FileName> &, const Options &o = Options());
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif