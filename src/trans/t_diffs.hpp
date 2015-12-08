#ifndef T_DIFFS_HPP
#define T_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TDiffs : public Analyzer
    {
        enum Assembler
        {
            Cuffdiffs,
            DESeq2,
            EdgeR,
        };
        
        enum RNALevel
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Assembler soft;

            // Only valid for Cuffdiffs
            RNALevel level;
        };

        struct Stats : public LinearStats, public MappingStats
        {
            Limit ss;

            // The keys depend on whether it's a gene or isoform analysis
            std::map<std::string, Counts> h;
        };

        static Stats analyze(const FileName &, const Options &o);
        static Stats analyze(const std::vector<DiffTest> &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif