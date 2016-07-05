#ifndef R_SAMPLE_HPP
#define R_SAMPLE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct RSample
    {
        enum class Method
        {
            Prop
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            // Defined if the method is proportion
            Proportion p;
            
            Method meth = Method::Prop;
        };

        struct Stats : public MappingStats
        {
            typedef std::vector<double> ChrTData;
            typedef std::vector<double> GenoData;
            
            // Statistics for the sequins
            ChrTData syn;
            
            // Statistics for the genome
            GenoData geno;
            
            /*
             * Coverage before subsampling
             */
            
            Coverage gBefore;
            Coverage sBefore;
            
            /*
             * Coverage after subsampling
             */
            
            Coverage gAfter;
            Coverage sAfter;
            
            // Proportion required for subsampling
            Proportion p = NAN;
        };

        static Stats stats(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif