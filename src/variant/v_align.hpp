#ifndef GI_V_ALIGN_HPP
#define GI_V_ALIGN_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct VAlign
    {
        typedef SingleMixtureOptions Options;
        
        struct Stats : public AlignmentStats
        {
            Performance p;

            // Distribution for the sequins
            SequinHist h = Standard::instance().r_var.hist();
        };

        static Stats report(const std::string &, const Options &o = Options());
    };
}

#endif