/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_DIFF_HPP
#define T_DIFF_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    class CountTable;
    
    struct TDiff : public Analyzer
    {
        template <typename Stats, typename Options> static Scripts writeCSV(const Stats &stats, const Options &o)
        {
            /*
             * Generating a file for differential results. The file should list the relevant information for
             * plotting an MA and LODR plot.
             *
             * We will need the following information:
             *
             *     - BaseMean
             *     - Expected LF
             *     - LogFold
             *     - LogFold SE
             *     - PValue
             */
            
            std::stringstream ss;
            ss << ",baseMean,elfc,lfc,lfcSE,pval";
            
//            for (const auto &i : stats.data)
            {
                const auto &ps        = stats.data.ps;
                const auto &ids       = stats.data.ids;
                const auto &logFs     = stats.data.logFs;
                const auto &eLogFs    = stats.data.eLogFs;
                const auto &logFSEs   = stats.data.logFSEs;
                const auto &baseMeans = stats.data.baseMeans;
                
                for (auto j = 0; j < ids.size(); j++)
                {
                    if (isnan(ps[j]))
                    {
                        ss << (boost::format("%1%,NA,NA,NA") % ids[j]).str();
                    }
                    else
                    {
                        ss << ((boost::format("%1%,%2%,%3%,%4%,%5%,%6%") % ids[j]
                                % toNA(baseMeans[j])
                                % toNA(eLogFs[j])
                                % toNA(logFs[j])
                                % toNA(logFSEs[j])
                                % toNA(ps[j])).str());
                    }
                }
            }
            
            return ss.str();
        }
        
        enum class Counting
        {
            None,
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

            Software dSoft;
            
            // Optional. Required for MA plot.
            Counting cSoft = Counting::None;
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
            
            Data data;

            // Average counts for each condition if provided
            std::vector<std::map<std::string, Counts>> avgs;
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
