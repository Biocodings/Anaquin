#ifndef V_KREPORT_HPP
#define V_KREPORT_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/v_allele.hpp"

namespace Anaquin
{
    struct VKReport : public Analyzer
    {
        struct Options : public AnalyzerOptions
        {
            FileName index;
        };

        struct Stats
        {
            VAllele::Stats allele;

            // The path for the analysis files
            Path output;            
        };
        
        static Stats analyze(const FileName &, const Options &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif