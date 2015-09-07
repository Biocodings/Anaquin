#ifndef GI_F_DISCOVER_HPP
#define GI_F_DISCOVER_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        struct Options : public SingleMixtureOptions
        {
            Options(Software soft = Software::TopHat, double fuzzy = 0) : soft(soft), fuzzy(fuzzy)
            {
                if (soft != TopHat && soft != Star)
                {
                    throw std::runtime_error("Only Tophat-Fusion and Star are supported");
                }
            }
            
            unsigned fuzzy;

            Software soft;
        };

        struct Stats : public LinearStats
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

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif