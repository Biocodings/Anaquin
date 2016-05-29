#ifndef M_KDIFF_HPP
#define M_KDIFF_HPP

#include "stats/analyzer.hpp"
#include "MetaQuin/MetaQuin.hpp"

namespace Anaquin
{
    struct MKDiff
    {
        // Represent a differential result for a MetaQuin sequin
        struct SequinDiff
        {
            inline bool operator<(const SequinDiff &x)  const { return id < x.id;      }
            inline bool operator==(const SequinDiff &x) const { return id != x.id;     }
            inline bool operator!=(const SequinDiff &x) const { return !operator==(x); }
            
            SequinID id;
            
            // Expected coverage for mixture A
            Coverage e1;
            
            // Expected coverage for mixture B
            Coverage e2;
            
            // Measured coverage for mixture A
            Coverage m1;
            
            // Measured coverage for mixture B
            Coverage m2;
            
            // Expected fold-change
            Coverage eFold;

            // Measured fold-change
            Coverage mFold;
        };

        struct Options : public DoubleMixtureOptions
        {
            // Required by the GCC compiler...
            Options() {}
            
            MAligner aligner;
            MSoftware soft;
        };

        struct Stats : public LinearStats, public MappingStats
        {
            // Absolute detection limit
            Limit limit;
            
            // Distribution of the sequins
            SequinHist hist;
            
            std::set<SequinID> miss;
            
            // Mapping for p-values
            std::map<SequinID, Probability> s2p;
            
            // Mapping for adjusted p-values
            std::map<SequinID, Probability> s2q;
            
            // Differential differences
            //std::set<SequinDiff> diffs;
        };

        static Stats analyze(const std::vector<FileName> &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
