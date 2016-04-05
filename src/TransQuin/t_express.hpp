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
            Isoform
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            // Gene or isoform?
            Metrics metrs;

            // What software generated the files?
            Software soft;
        };

        struct Stats : public MappingStats, public SequinStats
        {
            LinearStats data;
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

        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
