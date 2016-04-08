#ifndef M_DIFF_HPP
#define M_DIFF_HPP

#include "MetaQuin/m_blat.hpp"
#include "MetaQuin/m_abund.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MDiff
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

        struct Stats : public LinearStats, public MappingStats
        {
            // Alignment for the two alignment files
            MBlat::Stats align_1, align_2;
            
            // Absolute detection limit
            Limit limit;
            
            // Distribution of the sequins
            SequinHist hist;
            
            // Differential differences
            std::set<SequinDiff> diffs;
        };

        struct Options : public DoubleMixtureOptions
        {
            // Required by the GCC compiler...
            Options() {}
            
            // How the measured coverage is computed
            MAbund::CoverageMethod coverage = MAbund::WendySmooth;
            
            // An optional PSL file for mixture A
            FileName pA;

            // An optional PSL file for mixture B
            FileName pB;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
