#ifndef R_KEXPRESS_HPP
#define R_KEXPRESS_HPP

#include "stats/analyzer.hpp"
#include "RnaQuin/r_express.hpp"

namespace Anaquin
{
    struct RKReport : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        typedef AnalyzerOptions Options;
        
        struct Stats
        {
            Stats(const RExpress::Stats &stats) : stats(stats) {}
            
            RExpress::Stats stats;
            
            // Kallisto's abundance.tsv
            FileName abund;
        };
        
        static Stats analyze(const FileName &, const FileName &, const FileName &, const Options &);
        static void report(const FileName &, const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif