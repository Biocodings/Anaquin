#ifndef T_REPLICATE_HPP
#define T_REPLICATE_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        
        enum class Metrics
        {
            Gene,
            Isoform
        };
        
        enum class Software
        {
            Cufflinks,
            StringTie,
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

        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            if (files.size() < 2)
            {
                throw std::runtime_error("Two or more replicates required");
            }
            
            std::vector<TExpress::Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;
        }
        
        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);

        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
