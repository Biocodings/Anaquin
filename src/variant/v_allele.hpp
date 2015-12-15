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
            struct ChrT : public LinearStats, public MappingStats
            {
                Limit ss;
                
                Confusion m;
                
                long detected;
                
                // Sensitivity
                double sn;
                
                SequinHist h = Standard::instance().r_var.hist();
            };

            std::shared_ptr<ChrT> chrT;
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif