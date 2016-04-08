#ifndef L_EXPRESS_HPP
#define L_EXPRESS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct LExpress
    {
        typedef SingleMixtureOption Options;

        struct Stats : public AlignmentStats
        {
            struct Data : public LinearStats
            {
                // Histogram expected
                std::map<SequinID, Coverage> expect;
                
                // Histogram before normalization and adjustment
                SequinHist measured;
                
                // Histogram after normalization but before adjustment
                std::map<SequinID, Coverage> normalized;
                
                // Histogram after adjustment
                std::map<SequinID, Coverage> adjusted;
                
                // Adjusted abundance at the joined level
                std::map<LadderRef::JoinID, Coverage> joinAdjusted;
                
                // Expected size of the library
                Counts expTotal = 0;
                
                // Measured size of the library
                Counts obsTotal = 0;
            };
            
            Data data;

            // Histogram at the unjoined level
            SequinHist hist = Standard::instance().r_lad.hist();
            
            // Sensitivity at the joined level
            Limit s_joined;
            
            // Histogram at the joined level
            LadderRef::JoinHist h_joined  = Standard::instance().r_lad.joinHist();

            // Absolute detection limit
            Limit absolute;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif