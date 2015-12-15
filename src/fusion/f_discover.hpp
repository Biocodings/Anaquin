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

        struct Options : public AnalyzerOptions
        {
            Options(FDiscover::Aligner soft = Aligner::TopHat, double fuzzy = 0) : soft(soft), fuzzy(fuzzy)
            {
                if (soft != Aligner::TopHat && soft != Aligner::Star)
                {
                    throw std::runtime_error("Only Tophat-Fusion and Star are supported");
                }
            }
            
            unsigned fuzzy;

            FDiscover::Aligner soft;
        };

        struct Stats
        {
            struct ChrT : public LinearStats, public FusionStats
            {
                // Overall performance
                Confusion m;
                
                // Fraction of reference fusion detected
                double covered;
                
                // Distribution of the sequins
                SequinHist h = Standard::instance().r_fus.hist();            
            };

            std::shared_ptr<ChrT> chrT;
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif