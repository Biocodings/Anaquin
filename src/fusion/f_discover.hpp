#ifndef F_DISCOVER_HPP
#define F_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        enum Aligner
        {
            Star,
            TopHat,
        };

        struct Options : public FuzzyOptions
        {
            FDiscover::Aligner aligner;
        };

        struct Stats
        {
            struct Data : public LinearStats, public FusionStats
            {
                // Overall performance
                Confusion m;
                
                // Fraction of reference fusion detected
                double covered;
                
                // Distribution of the sequins
                SequinHist h = Standard::instance().r_fus.hist();            
            };

            std::map<ChromoID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif