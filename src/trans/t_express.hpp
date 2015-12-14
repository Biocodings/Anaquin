#ifndef T_ABUND_HPP
#define T_ABUND_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        enum Assembler
        {
            Cufflinks,
            StringTie,
        };
        
        enum RNALevel
        {
            Gene,
            Isoform
        };

        struct Stats : public LinearStats, public MappingStats
        {
            Limit ss;

            // The keys depend on whether it's a gene or isoform analysis
            std::map<std::string, Counts> h;
        };

        struct Options : public AnalyzerOptions
        {
            // This's required by gcc...
            Options() {}

            Assembler tool;
            
            RNALevel level = Isoform;
        };

        // Analyze for a single sample
        static Stats analyze(const FileName &, const Options &o);
        
        // Analyze for a single sample
        static Stats analyze(const std::vector<Expression> &, const Options &);

        // Analyze for multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
