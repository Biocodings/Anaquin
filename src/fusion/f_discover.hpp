#ifndef GI_F_DISCOVER_HPP
#define GI_F_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        enum Software
        {
            Star,
            TopHat,
        };

        struct Options : public AnalyzerOptions
        {
            Options(FDiscover::Software soft = Software::TopHat, double fuzzy = 0) : soft(soft), fuzzy(fuzzy)
            {
                if (soft != Software::TopHat && soft != Software::Star)
                {
                    throw std::runtime_error("Only Tophat-Fusion and Star are supported");
                }
            }
            
            unsigned fuzzy;

            FDiscover::Software soft;
        };

        struct Stats : public LinearStats, public FusionStats
        {
            // Overall performance
            Confusion m;

            // Fraction of reference fusion detected
            double covered;

            // Distribution of the sequins
            SequinHist h = Standard::instance().r_fus.hist();

            // Sequins failed to detect in the experiment
            MissingSequins miss;
        };

        static Stats report(const FileName &, const Options &o = Options());
    };
}

#endif