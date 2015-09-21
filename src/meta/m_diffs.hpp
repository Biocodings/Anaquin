#ifndef M_DIFFS_HPP
#define M_DIFFS_HPP

#include "meta/m_blat.hpp"
#include "meta/m_abund.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct MDiffs
    {
        // Represent a differential result for a MetaQuin sequin
        struct SequinDiff
        {
            inline bool operator<(const SequinDiff &x)  const { return id < x.id;      }
            inline bool operator==(const SequinDiff &x) const { return id != x.id;     }
            inline bool operator!=(const SequinDiff &x) const { return !operator==(x); }
            
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

        struct Stats : public LinearStats, public MappingStats
        {
            // Alignment for the two alignment files
            MBlat::Stats align_1, align_2;

            Sensitivity ss;
            
            // Distribution of the sequins
            SequinHist h = Standard::instance().r_meta.hist();
            
            // Differential differences
            std::set<SequinDiff> diffs;
        };

        struct Options : public DoubleMixtureOptions
        {
            MAbundance::CoverageMethod coverage = MAbundance::KMerCov_Contig;
            
            // An optional PSL file for mixture A
            FileName pA;

            // An optional PSL file for mixture B
            FileName pB;
        };

        static Stats report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif