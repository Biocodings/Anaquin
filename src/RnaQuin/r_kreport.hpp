#ifndef R_KREPORT_HPP
#define R_KREPORT_HPP

#include "stats/analyzer.hpp"
#include "RnaQuin/r_fold.hpp"
#include "RnaQuin/r_express.hpp"
#include "parsers/parser_exp.hpp"

namespace Anaquin
{
    struct RKReport : public Analyzer
    {
        struct Options : public AnalyzerOptions
        {
            FileName index;
        };

        struct Stats
        {
            // Kallsito quantified files
            std::vector<FileName> abunds;

            // Statistics for expression at the isoform level
            std::vector<RExpress::Stats> iExpress;
            
            // Statistics for expression at the gene level
            std::vector<RExpress::Stats> gExpress;
            
            // Statistics for differential at the isoform level
            RFold::Stats iFold;
            
            // Statistics for differential at the gene level
            RFold::Stats gFold;
            
            // The path for the analysis files
            Path output;

            ParserExp::Experiment exp;
        };
        
        static Stats analyze(const FileName &, const Options &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif