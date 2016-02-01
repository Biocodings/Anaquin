#ifndef V_ALLELE_HPP
#define V_ALLELE_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAllele
    {
        typedef FuzzyOptions Options;

        struct Stats
        {
            struct Data : public LinearStats, public MappingStats
            {
                Confusion m;
                
                long detected;
                
                // Sensitivity
                double sn;

                //SequinHist h = Standard::instance().r_var.hist();
            };

            Data data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif