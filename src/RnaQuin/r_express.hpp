/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_EXPRESS_HPP
#define R_EXPRESS_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMinor();

// Defined in resources.cpp
extern Anaquin::FileName MixRef();

// Defined in resources.cpp
extern Anaquin::FileName GTFRef();

namespace Anaquin
{
    struct TExpress : public Analyzer
    {
        /*
         * Generate summary statistics for a single sample
         */
        
        template <typename Stats, typename Options> static Scripts singleSummary(const Stats &stats,
                                                                                 const FileName &file,
                                                                                 const Units &units,
                                                                                 const Options &o)
        {
            return StatsWriter::inflectSummary(o.rAnnot,
                                               o.rAnnot,
                                               std::vector<FileName>     { file  },
                                               std::vector<SequinHist>   { stats.hist },
                                               std::vector<MappingStats> { stats },
                                               std::vector<LinearStats>  { stats },
                                               units);
        }
        
        /*
         * Generate summary statistics for multiple samples
         */
        
        template <typename Stats, typename Options> static Scripts multipleSummary(const std::vector<FileName> &files,
                                                                                   const std::vector<Stats>    &stats,
                                                                                   const Units &units,
                                                                                   const Options &o)
        {
            std::vector<SequinHist>   sHists;
            std::vector<LinearStats>  lStats;
            std::vector<MappingStats> mStats;
            
            for (auto i = 0; i < files.size(); i++)
            {
                mStats.push_back(stats[i]);
                sHists.push_back(stats[i].hist);
                lStats.push_back(stats[i]);
            }

            return StatsWriter::inflectSummary(o.rAnnot, o.rAnnot, files, sHists, mStats, lStats, units);
        }
        
        template <typename Stats> static Scripts multipleCSV(const std::vector<Stats> &stats)
        {
            std::set<SequinID> seqs;
            
            // This is the data structure that will be convenient
            std::map<unsigned, std::map<SequinID, Concent>> data;
            
            // Expected concentration
            std::map<SequinID, Concent> expect;
            
            std::stringstream ss;
            ss << "Seq\tExpected";
            
            for (auto i = 0; i < stats.size(); i++)
            {
                ss << ((boost::format("\tObserved%1%") % (i+1)).str());
                
                for (const auto &j : stats[i])
                {
                    seqs.insert(j.first);
                    expect[j.first]  = j.second.x;
                    data[i][j.first] = j.second.y;
                }
            }
            
            ss << "\n";
            
            for (const auto &seq : seqs)
            {
                ss << ((boost::format("%1%\t%2%") % seq % expect.at(seq)).str());
                
                for (auto i = 0; i < stats.size(); i++)
                {
                    if (data[i].count(seq))
                    {
                        ss << "\t" << data[i][seq];
                    }
                    else
                    {
                        ss << "\tNA";
                    }
                }
                
                ss << "\n";
            }
            
            return ss.str();
        }
        
        typedef ParserCufflink::Data TestData;
        
        enum class Metrics
        {
            Gene,
            Isoform
        };
        
        enum class Software
        {
            Kallisto,
            Cufflinks,
            StringTie,
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            // Gene or isoform?
            Metrics metrs;
        };
        
        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            // Empty Implementation
        };

        static Stats analyze(const FileName &, const Options &o);

        /*
         * 4. Generating major plot (but only if we have the isoforms...)
         */

        template <typename Options> static void generateRSplice(const FileName &output,
                                                                const FileName &csv,
                                                                const Options &o)
        {
            o.writer->open(output);
            o.writer->write(RWriter::createScript(csv, PlotTMinor()));
            o.writer->close();
        }

        /*
         * Generate an R-script for expected abundance against measured abundance
         */
        
        template <typename Options> static void generateRAbund(const FileName &output,
                                                               const FileName &csv,
                                                               const std::vector<Stats> &stats,
                                                               const Options &o)
        {
            o.info("Generating " + output);
            o.writer->open(output);
            
            const auto title = o.metrs == Metrics::Gene ? "Gene Expression" : "Isoform Expression";
            
            if (stats.size() == 1)
            {
                o.writer->write(RWriter::createScatterNeedLog(csv,
                                                              title,
                                                              "Expected Expression (log2)",
                                                              "Measured Expression (log2)",
                                                              "Expected",
                                                              "Measured", false));
            }
            else
            {
                o.writer->write(RWriter::createMultiScatterNeedLog(csv,
                                                                   title,
                                                                   "Expected Expression (log2)",
                                                                   "Measured Expression (log2)",
                                                                   "Expected",
                                                                   "Measured", false));
            }
            
            o.writer->close();
        }

        /*
         * Generate sequin statistics in the CSV format
         */
        
        template <typename Stats, typename Options> static void generateCSV(const FileName &output,
                                                                            const std::vector<Stats> &stats,
                                                                            const Options &o)
        {
            o.info("Generating " + output);
            o.writer->open(output);
            
            if (stats.size() == 1)
            {
                o.writer->write(StatsWriter::writeCSV(stats[0]));
            }
            else
            {
                o.writer->write(TExpress::multipleCSV(stats));
            }
            
            o.writer->close();
        }

        /*
         * Generate summary statistics for a single sample and multiple samples.
         */
        
        template <typename Stats, typename Options> static void generateSummary(const FileName &summary,
                                                                                const std::vector<FileName >&files,
                                                                                const std::vector<Stats> &stats,
                                                                                const Options  &o,
                                                                                const Units &units)
        {
            const auto &r = Standard::instance().r_trans;
            
            o.info("Generating " + summary);
            o.writer->open(summary);
            
            std::vector<SequinHist>   hists;
            std::vector<LinearStats>  lStats;
            std::vector<MappingStats> mStats;
            
            for (auto i = 0; i < files.size(); i++)
            {
                mStats.push_back(stats[i]);
                lStats.push_back(stats[i]);
                hists.push_back(stats[i].hist);
            }
            
            const auto title = o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed";

            const auto ms = StatsWriter::multiInfect(o.rAnnot, o.rAnnot, files, hists, mStats, lStats);

            const auto n_syn = toString(r.countGenesSyn()) + " " + units;
            const auto n_gen = toString(r.countGenesGen()) + " " + units;

            const auto format = "-------RnaExpression Output\n\n"
                                "Summary for input: %1%\n\n"
                                "       *Arithmetic average and standard deviation are shown\n\n"
                                "-------User Transcript Annotations\n\n"
                                "       Annotation file: %2%\n"
                                "       Synthetic: %3%\n"
                                "       Genome:    %4%\n\n"
                                "       Mixture file: %5%\n\n"
                                "-------%6%\n\n"
                                "       Synthetic: %7%\n"
                                "       Detection Sensitivity: %8% (attomol/ul) (%9%)\n\n"
                                "       Genome: %10%\n\n"
                                "-------Limit of Quantification (LOQ)\n"
                                "       *Estimated by piecewise segmented regression\n\n"
                                "       Break: %11% (%12%)\n\n"
                                "       *Below LOQ\n"
                                "       Intercept:   %13%\n"
                                "       Slope:       %14%\n"
                                "       Correlation: %15%\n"
                                "       R2:          %16%\n\n"
                                "       *Above LOQ\n"
                                "       Intercept:   %17%\n"
                                "       Slope:       %18%\n"
                                "       Correlation: %19%\n"
                                "       R2:          %20%\n\n"
                                "-------Linear regression (log2 scale)\n\n"
                                "       Correlation: %21%\n"
                                "       Slope:       %22%\n"
                                "       R2:          %23%\n"
                                "       F-statistic: %24%\n"
                                "       P-value:     %25%\n"
                                "       SSM:         %26%, DF: %27%\n"
                                "       SSE:         %28%, DF: %29%\n"
                                "       SST:         %30%, DF: %31%\n";

            o.writer->write((boost::format(format) % STRING(ms.files)
                                                   % GTFRef()
                                                   % n_syn    // 3
                                                   % n_gen    // 4
                                                   % MixRef() // 5
                                                   % title    // 6
                                                   % STRING(ms.n_syn)
                                                   % "????"
                                                   % "????"
                                                   % STRING(ms.n_gen)
                                                   % STRING(ms.b) // 11
                                                   % STRING(ms.bID)
                                                   % STRING(ms.lInt)
                                                   % STRING(ms.lSl)
                                                   % "????"
                                                   % STRING(ms.lR2) // 16
                                                   % STRING(ms.rInt)
                                                   % STRING(ms.rSl)
                                                   % "????"
                                                   % STRING(ms.rR2)
                                                   % STRING(ms.wLog.r) // 21
                                                   % STRING(ms.wLog.sl)
                                                   % STRING(ms.wLog.R2)
                                                   % STRING(ms.wLog.F)
                                                   % STRING(ms.wLog.p)
                                                   % STRING(ms.wLog.SSM) // 26
                                                   % STRING(ms.wLog.SSM_D)
                                                   % STRING(ms.wLog.SSE)
                                                   % STRING(ms.wLog.SSE_D)
                                                   % STRING(ms.wLog.SST)
                                                   % STRING(ms.wLog.SST_D) // 31
                             ).str());
            o.writer->close();
            
//            if (stats.size() == 1)
//            {
//                o.writer->write(singleSummary(stats[0], files[0], units, o));
//            }
//            else
//            {
//                o.writer->write(multipleSummary(files, stats, units, o));
//            }
        }
        
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<TExpress::Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;
        }

        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);

        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
