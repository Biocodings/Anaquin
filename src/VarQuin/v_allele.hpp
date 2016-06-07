#ifndef V_ALLELE_HPP
#define V_ALLELE_HPP

#include "stats/analyzer.hpp"
#include "VarQuin/VarQuin.hpp"

namespace Anaquin
{
    struct VAllele
    {
        enum class Input
        {
            VCF,
            Text,
            //Kallisto
        };
        
        struct Options : public AnalyzerOptions
        {
            Input input;
        };

        struct Stats : public MappingStats, public SequinStats, public VariantStats
        {
            // Statistics for all variants
            LinearStats all;
            
            /*
             * Not available for Kallisto
             */
            
            // Statistics for SNPs
            LinearStats snp;
            
            // Statistics for indels
            LinearStats ind;
            
            // Abundance for reference allele
            std::map<SequinID, Counts> readR;
            
            // Abundance for variant allele
            std::map<SequinID, Counts> readV;            
        };

        static Stats analyze(const FileName &, const Options &o = Options());
        static void report(const FileName &, const Options &o = Options());
    };
}

#endif