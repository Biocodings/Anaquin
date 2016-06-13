#ifndef R_SAMPLE_HPP
#define R_SAMPLE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RSample
    {
        enum class Software
        {
            None,
            Cufflinks,
        };
        
        /*
         * Specififes how subsampling should be done.
         */
        
        enum class Method
        {
            _1,
            _5,
            _10,
            _15,
            _20,
            _50
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Software soft;
            
            Method meth = Method::_10;
        };

        struct Stats : public MappingStats
        {
            typedef std::vector<double> ChrTData;
            typedef std::vector<double> GenoData;
            
            // Statistics for the sequins
            ChrTData chrT;
            
            // Statistics for the genome
            GenoData geno;
            
            Coverage genoBefore;
            Coverage chrTBefore;
            Coverage genoAfter;
            Coverage chrTAfter;
            
            // Proportion required for subsampling
            Proportion prop;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif