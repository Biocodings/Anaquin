#ifndef V_DISCOVER_HPP
#define V_DISCOVER_HPP

#include <vector>
#include "stats/analyzer.hpp"
#include "VARQuin/VARQuin.hpp"

namespace Anaquin
{
    struct VDiscover
    {
        struct Options : public AnalyzerOptions
        {
            Caller caller;
        };
        
        struct Stats : public MappingStats
        {
            struct ChrTData
            {
                // Expected allele frequency
                double eAllFreq;
                
                // Expected fold-change
                double eFold;
                
                CalledVariant query;

                // The sequin where it occurs
                const Variant *seq;
            };
            
            struct EndoData
            {
                CalledVariant query;
            };

            struct ChrTStats
            {
                std::vector<ChrTData> fps, tps;

                // Performance metrics
                Confusion m, m_snp, m_ind;
            };

            typedef std::vector<EndoData> EndoStats;
            
            // Statistics for synthetic variants
            ChrTStats chrT;

            // Statistics for endogenous variants
            EndoStats endo;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif