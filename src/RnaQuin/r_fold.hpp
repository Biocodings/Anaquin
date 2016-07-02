/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_FOLD_HPP
#define R_FOLD_HPP

#include "data/dtest.hpp"
#include <boost/format.hpp>
#include "stats/analyzer.hpp"

// Defined in resources.cpp
extern Anaquin::Scripts PlotTLODR();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTROC();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMA();

// Defined in resources.cpp
extern Anaquin::FileName MixRef();

// Defined in resources.cpp
extern Anaquin::FileName GTFRef();

namespace Anaquin
{
    struct RFold : public Analyzer
    {
        template <typename Stats, typename Options> static Scripts generateQuins(const Stats &stats, const Options &o)
        {
            std::stringstream ss;
            ss << "ID\tMean\tExpected\tMeasured\tSe\tPval\tQval\n";
            
            for (const auto &i : stats.data)
            {
                const auto &x = i.second;
                
                if (isnan(x.p))
                {
                    ss << (boost::format("%1%\tNA\tNA\tNA\tNA\n") % i.first).str();
                }
                else
                {
                    ss << ((boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\n") % i.first
                                                                                 % n2str(x.mean)
                                                                                 % n2str(x.exp)
                                                                                 % n2str(x.obs)
                                                                                 % n2str(x.se)
                                                                                 % p2str(x.q)
                                                                                 % p2str(x.p)).str());
                }
            }

            return ss.str();
        }

        template <typename Options> static void generateLODR(const FileName &file,
                                                             const FileName &csv,
                                                             const Options &o)
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTLODR()));
            o.writer->close();
        }

        template <typename Stats, typename Options> static void generateCSV(const FileName &file,
                                                                            const Stats &stats,
                                                                            const Options &o)
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write(RFold::generateQuins(stats, o));
            o.writer->close();
        }
        
        template <typename Stats, typename Options> static void generateSummary(const FileName &file,
                                                                                const FileName &src,
                                                                                const Stats &stats,
                                                                                const Options &o,
                                                                                const Units &units)
        {
            const auto &r = Standard::instance().r_trans;
            const auto lm = stats.linear(false);

            // No reference coordinate annotation given here
            const auto n_syn = o.metrs == Metrics::Gene ? r.countGeneSeqs() : r.countSeqs();

            const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");

            const auto summary = "-------RnaFoldChange Output\n\n"
                                 "       Summary for input: %1%\n\n"
                                 "-------Reference Annotations\n\n"
                                 "       Synthetic: %2% %3%\n"
                                 "       Mixture file: %4%\n\n"
                                 "-------%5%\n\n"
                                 "       Synthetic: %6% %3%\n"
                                 "       Genome:    %7% %3%\n\n"
                                 "       Detection Sensitivity: %8% (attomol/ul) (%9%)\n\n"
                                 "-------Linear regression (log2 scale)\n\n"
                                 "       Slope:       %10%\n"
                                 "       Correlation: %11%\n"
                                 "       R2:          %12%\n"
                                 "       F-statistic: %13%\n"
                                 "       P-value:     %14%\n"
                                 "       SSM:         %15%, DF: %16%\n"
                                 "       SSE:         %17%, DF: %18%\n"
                                 "       SST:         %19%, DF: %20%\n";
            o.generate(file);
            o.writer->open(file);
            o.writer->write((boost::format(summary) % src               // 1
                                                    % n_syn             // 2
                                                    % units             // 3
                                                    % MixRef()          // 4
                                                    % title             // 5
                                                    % stats.n_syn       // 6
                                                    % stats.n_gen       // 7
                                                    % stats.limit.abund // 8
                                                    % stats.limit.id    // 9
                                                    % lm.m              // 10
                                                    % lm.r              // 11
                                                    % lm.R2             // 12
                                                    % lm.F              // 13
                                                    % lm.p              // 14
                                                    % lm.SSM            // 15
                                                    % lm.SSM_D          // 16
                                                    % lm.SSE            // 17
                                                    % lm.SSE_D          // 18
                                                    % lm.SST            // 19
                                                    % lm.SST_D          // 20
                             ).str());
            o.writer->close();
        }
        
        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}
            Metrics metrs;
        };

        struct Stats : public LinearStats, public MappingStats, public AnalyzerStats
        {
            struct Data
            {
                // Expcted log-fold ratio
                Concent exp;
                
                // Measured log-fold ratio
                Concent obs;
                
                // Standard deviation
                double se;
                
                // Base mean
                double mean;
                
                Probability p, q;
            };
            
            std::map<SequinID, Data> data;
            
            /*
             * Optional inputs
             */
            
            // Normalized average counts for the replicates
            std::vector<double> means;
            
            // Log-fold ratios standard deviation
            std::vector<double> ses;
            
            // Average counts for each condition if provided
            //std::vector<std::map<std::string, Counts>> avgs;
        };

        /*
         * Classify and construct a vector of TP/FP/TN/FN, given the q-values and expected
         * fold-changes. The threshold for the TP classification is also required.
         */
        
        static std::vector<std::string> classify(const std::vector<double> &,
                                                 const std::vector<double> &,
                                                 double qCut,
                                                 double foldCut);

        static Stats analyze(const FileName &, const Options &o);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif