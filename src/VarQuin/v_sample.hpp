#ifndef V_SAMPLE_HPP
#define V_SAMPLE_HPP

#include "stats/analyzer.hpp"
#include "tools/coverage.hpp"

namespace Anaquin
{
    struct VSample
    {
        enum class Method
        {
            Mean,
            Median,
            ReadCount,
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
            
            Counts n_syn = 0;
            Counts n_gen = 0;
            
            /*
             * Fraction required to subsample in chrT. This works because chrT is a short
             * chromosome and almost certianly will have the highest coverage.
             */
            
            inline Proportion sample() const { return genC / synC; }
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}

            Method meth = Method::Mean;
        };

        static Stats stats(const FileName &file, const Options &o = Options());
        static void report(const FileName &file, const Options &o = Options());
    };
}

#endif