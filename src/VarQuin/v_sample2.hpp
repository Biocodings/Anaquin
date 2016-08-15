#ifndef V_SAMPLE2_HPP
#define V_SAMPLE2_HPP

#include "stats/analyzer.hpp"
#include "tools/coverage.hpp"

namespace Anaquin
{
//    struct SynGenAligns
//    {
//        Counts syn = 0;
//        Counts gen = 0;
//    };

    struct VSample2
    {
        enum class Method
        {
            Mean,
            Median,
            Reads,
            Prop,
        };
        
        struct Stats
        {
            // Intervals for synthetic
            ID2Intervals syn;
            
            // Intervals for genome
            ID2Intervals gen;
            
            // Raw coverage
            CoverageTool::Stats cov;
            
            // Calculated coverage for synthetic
            Coverage synC;
            
            // Calculated coverage for the query (eg: chr21)
            Coverage genC;
            
//            // Total alignments (not just sampling regions)
//            SynGenAligns tot;
//            
//            // Alignments within sampling regions
//            SynGenAligns samp;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Method meth = Method::Mean;
            
            // Defined only if meth==Prop
            Proportion p;
        };
        
        static Stats stats(const FileName &file, const Options &o = Options());
        static void report(const FileName &file, const Options &o = Options());
    };
}

#endif