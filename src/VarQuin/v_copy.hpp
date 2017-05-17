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
            
            CopyNumber copy = 2;
        };
        
        struct Stats
        {
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o);
        static void report  (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
