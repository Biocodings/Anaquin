#ifndef T_SAMPLE_HPP
#define T_SAMPLE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TSample
    {
        enum class Software
        {
            Cufflinks,
        };
        
        /*
         * Specififes how subsampling should be done.
         */
        
        enum class Method
        {
            // Calculate by averaging
            Mean,
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Software soft;
            
            Method meth = Method::Mean;
        };

        struct Stats : public MappingStats
        {
            typedef std::vector<double> ChrTData;
            typedef std::vector<double> GenoData;
            
            // Statistics for synthetic sequins
            ChrTData chrT;
            
            // Statistics for genomic transcripts
            GenoData geno;
            
            Depth genoBefore;
            Depth chrTBefore;
            Depth genoAfter;
            Depth chrTAfter;
            
            // Proportion required for subsampling
            Proportion prop;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif