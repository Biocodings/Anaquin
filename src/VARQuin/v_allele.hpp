#ifndef V_ALLELE_HPP
#define V_ALLELE_HPP

#include "stats/analyzer.hpp"
#include "VARQuin/VARQuin.hpp"

namespace Anaquin
{
    struct VAllele
    {
        struct Options : public AnalyzerOptions
        {
            Caller caller;
        };
        
        struct Stats : public MappingStats
        {
            typedef LinearStats ChrTData;

            struct ChrTStats
            {
                LinearStats tot, snp, ind;
            };
            
            typedef CalledVariant EndoData;
            typedef std::vector<EndoData> EndoStats;

            ChrTStats chrT;
            EndoStats endo;

            // Absolute detection limits
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif