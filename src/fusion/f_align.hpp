#ifndef GI_F_ALIGN_HPP
#define GI_F_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    class FAlign : public Analyzer
    {
        public:

            struct Stats
            {
            };
        
            struct Options : public SingleMixtureOptions
            {
                // Empty Implementation
            };

            static Stats report(const std::string &, const Options &options = Options());
    };
}

#endif