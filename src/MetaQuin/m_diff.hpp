/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef M_DIFF_HPP
#define M_DIFF_HPP

#include "stats/analyzer.hpp"
#include "MetaQuin/m_abund.hpp"

namespace Anaquin
{
    struct MDiff
    {
//        // Represent a differential result for a MetaQuin sequin
//        struct SequinDiff
//        {
//            inline bool operator<(const SequinDiff &x)  const { return id < x.id;      }
//            inline bool operator==(const SequinDiff &x) const { return id != x.id;     }
//            inline bool operator!=(const SequinDiff &x) const { return !operator==(x); }
//            
//            SequinID id;
//            
//            // Expected coverage for mixture A
//            Coverage e1;
//            
//            // Expected coverage for mixture B
//            Coverage e2;
//            
//            // Measured coverage for mixture A
//            Coverage m1;
//            
//            // Measured coverage for mixture B
//            Coverage m2;
//            
//            // Expected fold-change
//            Coverage eFold;
//            
//            // Measured fold-change
//            Coverage mFold;
//        };
        
        struct Options : public DoubleMixtureOptions
        {
            // Empty Implementation
        };
        
        struct Stats : public LinearStats, public MappingStats
        {
            MAbund::Stats stats1, stats2;
        };
        
        static Scripts generateRLinear(const FileName &, const Stats &, const Options &);
        
        static Stats analyze(const FileName &, const FileName &, const Options &o = Options());
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
