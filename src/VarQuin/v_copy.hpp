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
            // Estimated genomic normalization
            Proportion gNorm;
            
            double afterSeqs;
            
            VSample::GenomeSequins tBefore, tAfter;
            VSample::GenomeSequins sBefore, sAfter;
            
            VSample::CalibrateStats before;
            ParserBAMBED::Stats after;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o);
        static void report  (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
