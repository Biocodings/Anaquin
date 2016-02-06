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
            
            // Significance level
            Probability sign = 0.1;
        };

        struct Stats : public MappingStats
        {
            typedef VariantMatch ChrTData;
            
            struct ChrTStats
            {
                std::vector<ChrTData> fps, tps, tns, fns;

                // Performance metrics
                Confusion m, m_snp, m_ind;
            };

            typedef CalledVariant EndoData;
            
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