#ifndef R_KEXPRESS_HPP
#define R_KEXPRESS_HPP

#include "stats/analyzer.hpp"
#include "RnaQuin/r_fold.hpp"
#include "RnaQuin/r_express.hpp"

namespace Anaquin
{
    struct RKReport : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        
        struct Options : public AnalyzerOptions
        {
            FileName index;
        };

        struct Stats
        {
            // Kallsito quantification files
            std::vector<FileName> kFiles;
            
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
        };
        
        static Stats analyze(const std::vector<FileName> &, const Options &);        
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif