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

        struct Stats
        {
            // Overall performance
            Confusion m;

            // Fraction of reference fusion detected
            double covered;

            // Linear model of the coverage
            LinearStats cov;

            // Distribution of the sequins
            SequinHist h;// TODO: = Analyzer::seqHist();

            // Sequins failed to detect in the experiment
            MissingSequins miss;
        };

        static Stats analyze(const std::string &, const Options &options = Options());
    };
}

#endif