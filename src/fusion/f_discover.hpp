#ifndef F_DISCOVER_HPP
#define F_DISCOVER_HPP

#include "fusion/FUSQuin.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct FDiscover
    {
        struct Options : public FuzzyOptions
        {
            FusionCaller caller;
        };

        struct Stats : public FusionStats, public SequinStats
        {
            struct Data
            {
                inline Counts countFP()    const { return fps.size(); }
                inline Counts countKnown() const { return known;      }
                
                // Proportion of fusions detected
                //Proportion covered;

                // Number of detected fusions
                Counts detect;

                // Number of known fusions
                Counts known;

                // Distribution of the fusion
                SequinHist hist;
                
                // List of true-positive fusions
                std::vector<const FusionRef::KnownFusion *> tps;
                
                // List of false-positive fusions
                std::vector<FUSQuin::FalsePositive> fps;
            };

            std::map<ChromoID, Data> data;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif