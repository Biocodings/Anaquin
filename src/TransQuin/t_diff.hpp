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
        
        enum class Software
        {
            Cuffdiff,
            edgeR,
            DESeq2,
        };
        
        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}

            // Files for the count table
            std::vector<FileName> counts;
            
            Metrics metrs = Metrics::Gene;

            Software  dSoft;
            CountSoft cSoft;
        };

        struct Stats : public MappingStats, public SequinStats
        {
            struct Data : public LinearStats
            {
                // Detected features
                std::vector<FeatureID> ids;

                // Probability under the null hypothesis
                std::vector<double> ps;
                
                // Log-fold ratios
                std::vector<double> logFs;
                
                // Expected log-fold ratios (only for the synthetic)
                std::vector<double> eLogFs;
                
                /*
                 * Optional inputs. For example, Cuffdiffs wouldn't give them.
                 */

                // Normalized average counts for the replicates
                std::vector<double> baseMeans;

                // Log-fold ratios standard deviation
                std::vector<double> logFSEs;
            };
            
            std::map<ChrID, Data> data;

            // Count table
            std::shared_ptr<CountTable> counts;

            // Average counts for each condition if provided
            std::vector<std::map<std::string, Counts>> avgs;

            // Detection limit
            Limit limit;
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
