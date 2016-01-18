/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_DIFFS_HPP
#define T_DIFFS_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    class CountTable;
    
    struct TDiffs : public Analyzer
    {
        enum class CountSoft
        {
            HTSeqCount,
        };
        
        enum class DiffSoft
        {
            Cuffdiff,
            edgeR,
            DESeq2,
        };
        
        enum class Level
        {
            Gene,
            Exon,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}

            // Files for the count table
            std::vector<FileName> counts;
            
            // Default to gene level
            Level lvl = Level::Gene;

            DiffSoft  dSoft;
            CountSoft cSoft;
        };

        struct Stats : public MappingStats
        {
            struct Data : public LinearStats
            {
                // Detected features
                std::vector<GenericID> ids;

                // Probability under the null hypothesis
                std::vector<double> ps;
                
                // Q-probability (controlled for multiple testing)
                std::vector<double> qs;
                
                // Log-fold changes
                std::vector<double> logFs;
            };
            
            std::map<ChromoID, Data> data;

            // Average counts for each condition if provided
            std::vector<std::map<std::string, Counts>> avgs;
            
            std::shared_ptr<CountTable> counts;
            
            // Detection limit            
            Limit limit;

            inline bool haveCounts() const { return !avgs.empty(); }
        };

        /*
         * Classify and construct a vector of TP/FP/TN/FN, given a vector of q-values and expected fold-changes. The threshold
         * for the TP classification is also required.
         */
        
        static std::vector<std::string> classify(const std::vector<double> &, const std::vector<double> &, double qCut, double foldCut);
        
        static Stats analyze(const FileName &, const Options &o);
        static Stats analyze(const std::vector<DiffTest> &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif
