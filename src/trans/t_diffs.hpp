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
            Assembler soft = Assembler::Cuffdiffs;

            // Only valid for Cuffdiffs
            RNALevel level;
        };

        struct Stats
        {
            struct ChrT : public LinearStats, public MappingStats
            {
                Limit ss;
                
                // The keys depend on whether it's a gene or isoform analysis
                std::map<std::string, Counts> h;
            };
            
            struct Gencode : public LinearStats, public MappingStats
            {
                
            };
            
            std::shared_ptr<ChrT> chrT;
            std::shared_ptr<Gencode> gcode;
        };

        // Analyze a single sample
        static Stats analyze(const FileName &, const Options &o);
        
        // Analyze a single sample
        static Stats analyze(const std::vector<DiffTest> &, const Options &o);

        // Analyze multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif