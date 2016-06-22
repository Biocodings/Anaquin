#ifndef V_FREQUENCY_HPP
#define V_FREQUENCY_HPP

#include "stats/analyzer.hpp"
#include "tools/vcf_data.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VFreq
    {
        typedef VarInput Input;

        struct Options : public AnalyzerOptions
        {
            Input input;
        };

        struct Stats : public MappingStats
        {
            struct Data
            {
                // Measured minor allele frequency
                Proportion af;
            };
            
            std::map<ChrID, Data> data;
            
            VCFData vData;

            // Statistics for all variants
            LinearStats vars;

            // Statistics for SNPs
            LinearStats snp;
            
            // Statistics for indels
            LinearStats ind;
            
            std::map<long, Counts> readR;
            std::map<long, Counts> readV;
            std::map<long, Counts> depth;
            
            // Distribution for the variants
            std::map<ChrID, HashHist> hist;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif