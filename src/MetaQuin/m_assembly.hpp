#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MAssembly
    {
        struct Stats : public LinearStats, public DAsssembly::DenoAssemblyImpl
        {
            inline SequinID findSeq(const ContigID &cID) const
            {
                auto x = cID;
                
                if (soft == MSoftware::RayMeta && aligner == MAligner::Blat)
                {
                    /*
                     * BLAT doesn't do spacing. For example, "contig-0 52976 nucleotides"
                     * would be given as "contig-0".
                     */
                    
                    x = Tokens::first(cID, " ");
                }
                
                return c2s.count(x) ? c2s.at(x) : "";
            }
            
            inline Proportion covered() const
            {
                return static_cast<Proportion>(match) / (match + mismatch);
            }
            
            /*
             * Usually, they're not needed but we need them here to check for special cases.
             */
            
            MAligner aligner;
            MSoftware soft;
            
            // Total mismatching bases
            Base mismatch = 0;
            
            // Total matching bases
            Base match = 0;
            
            std::map<SequinID, std::vector<ContigID>> s2c;
            
            // Mapping from contigs to sequins
            std::map<ContigID, SequinID> c2s;
            
            // Mapping from contigs to their length
            std::map<ContigID, Base> c2l;
            
            // Mapping from contigs to number of bases assembled
            std::map<ContigID, Base> c2a;
            
            // Statistics for de-novo assembly
            DAsssembly::Stats<DAsssembly::Contig> dnovo;
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            MAligner aligner;
            MSoftware soft;
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &o = Options());
        static void  report (const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
