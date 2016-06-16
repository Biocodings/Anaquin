/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_DIFF_HPP
#define R_DIFF_HPP

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
    class CountTable;
    
    struct RDiff : public Analyzer
    {
        template <typename Stats, typename Options> static Scripts writeCSV(const Stats &stats, const Options &o)
        {
            /*
             * Generating a file for differential analysis. The file should give the relevant data
             * for MA and LODR plot.
             */
            
            std::stringstream ss;
            ss << "Seq\tMean\tExpected\tMeasured\tSe\tPval\tQval\n";
            
            const auto &ps    = stats.ps;
            const auto &qs    = stats.qs;
            const auto &ids   = stats.ids;
            const auto &elfs  = stats.elfs;
            const auto &mlfs  = stats.mlfs;
            const auto &ses   = stats.ses;
            const auto &means = stats.means;
            
            for (auto j = 0; j < ids.size(); j++)
            {
                if (Standard::isSynthetic(stats.cIDs[j]))
                {
                    if (isnan(ps[j]))
                    {
                        ss << (boost::format("%1%\tNA\tNA\tNA\tNA\n") % ids[j]).str();
                    }
                    else
                    {
                        ss << ((boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\n") % ids[j]
                                                                                     % n2str(means[j])
                                                                                     % n2str(elfs[j])
                                                                                     % n2str(mlfs[j])
                                                                                     % n2str(ses[j])
                                                                                     % p2str(ps[j])
                                                                                     % p2str(qs[j])).str());
                    }
                }
            }
            
            return ss.str();
        }
        
        template <typename Options> static void generateMA(const FileName &file,
                                                           const FileName &csv,
                                                           const Options &o)
        {
            o.generate(file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTMA()));
            o.writer->close();
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
            o.writer->write(RDiff::writeCSV(stats, o));
            o.writer->close();
        }
        
        template <typename Stats, typename Options> static void generateSummary(const FileName &file,
                                                                                const Stats &stats,
                                                                                const Options &o,
                                                                                const Units &units)
        {
            const auto &r = Standard::instance().r_trans;

            const auto lm = stats.linear(true);
            const auto summary = "-------RnaFoldChange Output\n\n"
                                 "       Reference mixture file: %1%\n"
                                 "       User fold-change file:  %2%\n\n"
                                 "-------User Transcript Annotations\n\n"
                                 "       Annotation file: %3%\n"
                                 "       Synthetic: %4%\n"
                                 "       Genome:    %5%\n\n"
                                 "-------Genes Expressed\n\n"
                                 "       Synthetic: %6%\n"
                                 "       Detection Sensitivity: %7% (attomol/ul) (%8%)\n\n"
                                 "       Genome:    %9%\n\n"
                                 "-------Linear regression (log2 scale)\n\n"
                                 "       Correlation: %10%\n"
                                 "       Slope:       %11%\n"
                                 "       R2:          %12%\n"
                                 "       F-statistic: %13%\n"
                                 "       P-value:     %14%\n"
                                 "       SSM:         %15%, DF: %16%\n"
                                 "       SSE:         %17%, DF: %18%\n"
                                 "       SST:         %19%, DF: %20%\n";
            o.generate(file);
            o.writer->open(file);
            o.writer->write((boost::format(summary) % MixRef()
                                                    % file
                                                    % GTFRef()
                                                    % r.countGeneSyn()
                                                    % r.countGeneGen()
                                                    % stats.n_syn
                                                    % stats.limit.abund
                                                    % stats.limit.id
                                                    % stats.n_gen
                                                    % lm.r
                                                    % lm.m
                                                    % lm.R2
                                                    % lm.F
                                                    % lm.p
                                                    % lm.SSM
                                                    % lm.SSM_D
                                                    % lm.SSE
                                                    % lm.SSE_D
                                                    % lm.SST
                                                    % lm.SST_D).str());
            o.writer->close();
        }

        enum class Counting
        {
            None,
            HTSeqCount,
        };
        
        enum class Software
        {
            edgeR,
            DESeq2,
            Sleuth,
            Cuffdiff,
        };
        
        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}
            Metrics metrs = Metrics::Gene;
        };

        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            // Detected features (genes or isoforms)
            std::vector<FeatureID> ids;
            
            // Probability under the null hypothesis
            std::vector<Probability> ps, qs;

            // Expected log-fold ratios
            std::vector<Concent> elfs;

            std::vector<ChrID> cIDs;
            
            // Measured log-fold ratios
            std::vector<Concent> mlfs;

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
        static Stats analyze(const std::vector<DiffTest> &, const Options &o);

        static void report(const FileName &, const Options &o = Options());
    };
}

#endif