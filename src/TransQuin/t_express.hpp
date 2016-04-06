#ifndef T_EXPRESS_HPP
#define T_EXPRESS_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        
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
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
