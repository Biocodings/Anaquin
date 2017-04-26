#ifndef V_TRIM_HPP
#define V_TRIM_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VTrim
    {
        enum class Method
        {
            Left,
            Right,
            LeftRight
        };
        
        struct Stats
        {
            // Number of trimmed reads on left
            Counts left = 0;
            
            // Number of trimmed reads on right
            Counts right = 0;
            
            // Number of alignments before trimming
            Counts before = 0;
            
            // Number of alignments after trimming
            Counts after  = 0;
            
            // Number of reference regions
            Counts nRegs;
            
            // Reads trimmed on left and right
            std::set<ReadName> lTrim, rTrim;
        };

        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            // How to trim?
            Method meth = Method::LeftRight;
            
            // How much edge effects?
            Base trim = 1;
        };
        
        static Stats analyze(const FileName &, const Options &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
