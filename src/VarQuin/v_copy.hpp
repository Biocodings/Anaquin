#ifndef V_COPY_HPP
#define V_COPY_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/v_sample.hpp"

namespace Anaquin
{
    struct VCopy
    {
        typedef unsigned CopyNumber;
        
        struct Options : public VSample::Options
        {
            Options() {}
            
            // Copy number that we want to calibrate to the genome
            CopyNumber gen = 2;
        };

        struct Stats : public SequinStats, public LimitStats
        {
            // Estimated normalization (genome)
            Proportion gNorm;
            
            VSample::CalibrateStats before;
            ParserBAMBED::Stats after;
            
            // Average coverage for the genome (2n)
            Coverage gEndo;
            
            // Average coverage before subsampling (2n)
            Coverage bSeqs;
            
            // Average coverage after subsampling (2n)
            Coverage aSeqs;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o);
        static void report  (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
