#ifndef T_ASSEMBLY_HPP
#define T_ASSEMBLY_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TAssembly : Analyzer
    {
        struct Options : FuzzyOptions
        {
            // Path for the reference and query GTF
            FileName ref, query;
        };

        struct Stats
        {
            struct ChrT : public MappingStats
            {
                SequinHist hb = Standard::instance().r_trans.geneHist("chrT");
                SequinHist he = Standard::instance().r_trans.hist();
                SequinHist hi = Standard::instance().r_trans.hist();
                SequinHist ht = Standard::instance().r_trans.hist();
                
                double exonSP,   exonSN;
                double baseSP,   baseSN;
                double transSP,  transSN;
                double intronSP, intronSN;
                
                Limit sb, si, se, st;
            };
            
            struct Experiment : public MappingStats
            {
                
            };
            
            std::shared_ptr<ChrT> chrT;
            std::shared_ptr<Experiment> expT;
        };

        // Analyze a single sample
        static Stats analyze(const FileName &, const Options &o = Options());
        
        // Analyze multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &, const Options &o = Options());

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif