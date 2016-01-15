#ifndef T_EXPRESS_HPP
#define T_EXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        enum class Software
        {
            Cufflinks,
            StringTie,
        };
        
        enum class Metrics
        {
            Gene,
            Exon,
            Isoform
        };

        struct Stats : public MappingStats
        {
            typedef LinearStats Data;

            std::map<ChromoID, Data> data;
            
            Limit limit;
        };

        struct Options : public AnalyzerOptions
        {
            // This's required by gcc...
            Options() {}

            Metrics metrs = Metrics::Gene;
            
            Software soft;
        };

        static Stats analyze(const FileName &, const Options &o);
        static Stats analyze(const std::vector<Expression> &, const Options &);

        // Analyze for multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<TExpress::Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;
        }

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
