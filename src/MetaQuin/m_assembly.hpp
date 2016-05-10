#ifndef M_ASSEMBLY_HPP
#define M_ASSEMBLY_HPP

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

        struct Stats : public LinearStats, public DAsssembly::DenoAssemblyImpl
        {
            inline SequinID findSeq(const ContigID &cID) const
            {
                return c2s.count(cID) ? c2s.at(cID) : "";
            }
            
            std::map<ContigID, SequinID> c2s;
            std::map<SequinID, std::vector<ContigID>> s2c;
            
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