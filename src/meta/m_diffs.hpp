#ifndef GI_M_DIFFS_HPP
#define GI_M_DIFFS_HPP

#include "meta/m_blat.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MDiffs
    {
        // Represents a differential result for a metaquin
        struct SequinDiff
        {
            SequinID id;
            
            // Expected coverage for mixture A
            Coverage ex_A;
            
            // Expected coverage for mixture B
            Coverage ex_B;
            
            // Observed coverage for mixture A
            Coverage ob_A;
            
            // Observed coverage for mixture B
            Coverage ob_B;
            
            // Expected fold-change
            Coverage ex_fold;

            // Observed fold-change
            Coverage ob_fold;
        };

        struct Stats : public LinearStats
        {
            // Alignment for the two alignment files
            MBlat::Stats align_1, align_2;
            
            std::vector<SequinDiff> diffs;
        };

        struct Options : public DoubleMixtureOptions
        {
            // An optional PSL file for mixture A
            std::string pA;

            // An optional PSL file for mixture B
            std::string pB;
        };

        static Stats report(const std::string &, const std::string &, const Options &options = Options());
    };
}

#endif