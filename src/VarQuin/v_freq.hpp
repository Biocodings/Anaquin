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

        struct Stats : public MappingStats, public SequinStats
        {
            // Statistics for all variants
            LinearStats vars;
            
            struct Data
            {
                // Measured minor allele frequency
                Proportion af;
            };
            
            std::map<ChrID, Data> data;
            
            VCFData vData;

            /*
             * Not available for Kallisto
             */
            
            // Statistics for SNPs
            LinearStats snp;
            
            // Statistics for indels
            LinearStats ind;
            
            // Reads for reference allele
            std::map<SequinID, Counts> readR;
            
            // Reads for variant allele
            std::map<SequinID, Counts> readV;
            
            // Distribution for the variants
            std::map<ChrID, HashHist> hist;
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif