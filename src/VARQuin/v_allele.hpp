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
        
        struct Stats : public MappingStats, public SequinStats, public VariantStats
        {
            struct ChrTData
            {
                LinearStats tot, snp, ind;
            };

            typedef std::vector<CalledVariant> EndoData;

            ChrTData chrT;
            EndoData endo;

            // Absolute detection limit
            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif