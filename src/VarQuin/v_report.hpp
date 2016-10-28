#ifndef V_REPORT_HPP
#define V_REPORT_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/v_allele.hpp"
#include "parsers/parser_exp.hpp"

namespace Anaquin
{
    struct VReport : public Analyzer
    {
        struct Options : public AnalyzerOptions
        {
            FileName index;
        };

        struct Stats
        {
            // Allele frequency
            VAllele::Stats allele;

            ParserExp::Experiment exp;
        };
        
        static Stats analyze(const FileName &, const Options &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
