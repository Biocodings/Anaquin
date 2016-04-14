/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_DIFF_HPP
#define T_DIFF_HPP

#include "data/dtest.hpp"
#include "stats/analyzer.hpp"

// Defined in resources.cpp
extern Anaquin::Scripts PlotTFold();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTLODR();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTROC();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMA();

namespace Anaquin
{
    class CountTable;
    
    struct TDiff : public Analyzer
    {
        template <typename Stats, typename Options> static Scripts writeCSV(const Stats &stats, const Options &o)
        {
            /*
             * Generating a file for differential analysis. The file should give the relevant data
             * for MA and LODR plot.
             */
            
            std::stringstream ss;
            ss << "sequin\tmean\texpected\tmeasured\tse\tpval\tqval\n";
            
            const auto &ps        = stats.ps;
            const auto &qs        = stats.qs;
            const auto &ids       = stats.ids;
            const auto &logFs     = stats.mLogFs;
            const auto &eLogFs    = stats.eLogFs;
            const auto &logFSEs   = stats.logFSEs;
            const auto &baseMeans = stats.baseMeans;
            
            for (auto j = 0; j < ids.size(); j++)
            {
                if (isnan(ps[j]))
                {
                    ss << (boost::format("%1%\tNA\tNA\tNA\tNA\n") % ids[j]).str();
                }
                else
                {
                    ss << ((boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\n") % ids[j]
                                                                                 % toNA(baseMeans[j])
                                                                                 % toNA(eLogFs[j])
                                                                                 % toNA(logFs[j])
                                                                                 % toNA(logFSEs[j])
                                                                                 % toNA(ps[j])
                                                                                 % toNA(qs[j])).str());
                }
            }
            
            return ss.str();
        }
        
        template <typename Options> static void generateMA(const FileName &file,
                                                           const FileName &csv,
                                                           const Options &o)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTMA()));
            o.writer->close();
        }
        
        template <typename Options> static void generateLODR(const FileName &file,
                                                             const FileName &csv,
                                                             const Options &o)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTLODR()));
            o.writer->close();
        }

        template <typename Options> static void generateROC(const FileName &file,
                                                            const FileName &csv,
                                                            const Options &o)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTROC()));
            o.writer->close();
        }
        
        template <typename Options> static void generateFoldR(const FileName &file,
                                                              const FileName &csv,
                                                              const Options &o)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(RWriter::createScript(csv, PlotTFold()));
            o.writer->close();
        }

        template <typename Stats, typename Options> static void generateCSV(const FileName &file,
                                                                            const Stats &stats,
                                                                            const Options &o)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(TDiff::writeCSV(stats, o));
            o.writer->close();
        }
        
        template <typename Stats, typename Options> static void generateSummary(const FileName &file,
                                                                                const Stats &stats,
                                                                                const Options &o,
                                                                                const Units &units)
        {
            o.info("Generating " + file);
            o.writer->open(file);
            o.writer->write(StatsWriter::linearSummary(file, o.rChrT, stats, stats, stats.hist, units));
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

            // Files for the count table
            std::vector<FileName> counts;
            
            Metrics metrs = Metrics::Gene;

            Software dSoft;
            
            // Optional. Required for MA plot.
            Counting cSoft = Counting::None;
        };

        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            // Detected features (genes or isoforms)
            std::vector<FeatureID> ids;
            
            // Probability under the null hypothesis
            std::vector<Probability> ps, qs;
            
            // Expected log-fold ratios
            std::vector<Concent> eLogFs;

            // Measured log-fold ratios
            std::vector<Concent> mLogFs;

            /*
             * Optional inputs. For example, Cuffdiffs wouldn't give them.
             */
            
            // Normalized average counts for the replicates
            std::vector<double> baseMeans;
            
            // Log-fold ratios standard deviation
            std::vector<double> logFSEs;
            
            // Average counts for each condition if provided
            std::vector<std::map<std::string, Counts>> avgs;
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