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
    struct RExpress : public Analyzer
    {
        template <typename Stats> static Scripts multipleCSV(const std::vector<Stats> &stats)
        {
            std::set<SequinID> seqs;
            
            // This is the data structure that will be convenient
            std::map<unsigned, std::map<SequinID, Concent>> data;
            
            // Expected concentration
            std::map<SequinID, Concent> expect;
            
            std::stringstream ss;
            ss << "ID\tExpected";
            
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
        
        enum class Inputs
        {
            GTF,
            Text
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}
            
            Inputs inputs;

            // Gene or isoform?
            Metrics metrs;
        };
        
        struct GenData
        {
            // Eg: FPKM
            double abund = NAN;
        };
        
        struct Stats : public LinearStats, public MappingStats, public SequinStats
        {
            std::map<GenoID, GenData> gData;
        };

        static Stats analyze(const FileName &, const Options &o);

        /*
         * Generating major plot (but only if we have the isoforms...)
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
        
        template <typename Options> static void generateR(const FileName &output,
                                                          const FileName &csv,
                                                          const std::vector<Stats> &stats,
                                                          const Options &o)
        {
            o.info("Generating " + output);
            o.writer->open(output);
            
            const auto title = o.metrs == Metrics::Gene ? "Gene Expression" : "Isoform Expression";
            
            if (stats.size() == 1)
            {
                o.writer->write(RWriter::createScatterNeedLog(csv, title,
                                                              "Expected Expression (log2)",
                                                              "Measured Expression (log2)",
                                                              "Expected",
                                                              "Measured", true));
            }
            else
            {
                o.writer->write(RWriter::createMultiScatter(csv, title,
                                                            "Expected Expression (log2)",
                                                            "Measured Expression (log2)",
                                                            "Expected",
                                                            "Measured", true, true));
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
                o.writer->write(RExpress::multipleCSV(stats));
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
            
            // Detection limit for the replicates
            Limit limit;

            for (auto i = 0; i < files.size(); i++)
            {
                mStats.push_back(stats[i]);
                lStats.push_back(stats[i]);
                hists.push_back(stats[i].hist);

                if (isnan(limit.abund) || stats[i].limit.abund < limit.abund)
                {
                    limit = stats[i].limit;
                }
            }
            
            const auto title = (o.metrs == Metrics::Gene ? "Genes Expressed" : "Isoform Expressed");

            /*
             * TODO: Need to do it in the log scale!
             */
            
            const auto ms = StatsWriter::multiInfect(o.rAnnot, o.rAnnot, files, hists, mStats, lStats);

            const auto n_syn = toString(r.countGeneSyn()) + " " + units;
            const auto n_gen = toString(r.countGeneGen()) + " " + units;

            // Breakpoint estimated by piecewise regression
            const auto b = ms.b.mean();

            // Number of genomic features above the breakpoint
            SReals n_above;

            // Number of genomic features below the breakpoint
            SReals n_below;
            
            for (const auto &i : stats)
            {
                Counts above = 0;
                Counts below = 0;

                for (const auto &j : i.gData)
                {
                    assert(!isnan(j.second.abund));
                    
                    if (j.second.abund >= b)
                    {
                        above++;
                    }
                    else
                    {
                        below++;
                    }
                }
                
                n_above.add(above);
                n_below.add(below);
            }

            const auto format = "-------RnaExpression Output\n"
                                "       Summary for input: %1%\n"
                                "       *Arithmetic average and standard deviation are shown\n\n"
                                "-------User Transcript Annotations\n\n"
                                "       Synthetic: %2%\n"
                                "       Genome:    %3%\n\n"
                                "       Mixture file: %4%\n\n"
                                "-------%5%\n\n"
                                "       Synthetic: %6%\n"
                                "       Detection Sensitivity: %7% (attomol/ul) (%8%)\n\n"
                                "       Genome: %9%\n\n"
                                "-------Limit of Quantification (LOQ)\n"
                                "       *Estimated by piecewise segmented regression\n\n"
                                "       Break: %10% (%11%) attomol/ul\n\n"
                                "       *Below LOQ\n"
                                "       Intercept:   %12%\n"
                                "       Slope:       %13%\n"
                                "       Correlation: %14%\n"
                                "       R2:          %15%\n"
                                "       Genome:      %16%\n\n"
                                "       *Above LOQ\n"
                                "       Intercept:   %17%\n"
                                "       Slope:       %18%\n"
                                "       Correlation: %19%\n"
                                "       R2:          %20%\n"
                                "       Genome:      %21%\n\n"
                                "-------Linear regression (log2 scale)\n\n"
                                "       Correlation: %22%\n"
                                "       Slope:       %23%\n"
                                "       R2:          %24%\n"
                                "       F-statistic: %25%\n"
                                "       P-value:     %26%\n"
                                "       SSM:         %27%, DF: %28%\n"
                                "       SSE:         %29%, DF: %30%\n"
                                "       SST:         %31%, DF: %32%\n";

            o.writer->write((boost::format(format) % STRING(ms.files)      // 1
                                                   % n_syn                 // 2
                                                   % n_gen                 // 3
                                                   % MixRef()              // 4
                                                   % title                 // 5
                                                   % STRING(ms.n_syn)      // 6
                                                   % limit.abund           // 7
                                                   % limit.id              // 8
                                                   % STRING(ms.n_gen)      // 9
                                                   % STRING(ms.b)          // 10
                                                   % STRING(ms.bID)        // 11
                                                   % STRING(ms.lInt)       // 12
                                                   % STRING(ms.lSl)        // 13
                                                   % STRING(ms.lr)         // 14
                                                   % STRING(ms.lR2)        // 15
                                                   % STRING(n_below)       // 16
                                                   % STRING(ms.rInt)       // 17
                                                   % STRING(ms.rSl)        // 18
                                                   % STRING(ms.rr)         // 19
                                                   % STRING(ms.rR2)        // 20
                                                   % STRING(n_above)       // 21
                                                   % STRING(ms.wLog.r)     // 22
                                                   % STRING(ms.wLog.sl)    // 23
                                                   % STRING(ms.wLog.R2)    // 24
                                                   % STRING(ms.wLog.F)     // 25
                                                   % STRING(ms.wLog.p)     // 26
                                                   % STRING(ms.wLog.SSM)   // 27
                                                   % STRING(ms.wLog.SSM_D) // 28
                                                   % STRING(ms.wLog.SSE)   // 29
                                                   % STRING(ms.wLog.SSE_D) // 30
                                                   % STRING(ms.wLog.SST)   // 31
                                                   % STRING(ms.wLog.SST_D) // 32
                             ).str());
            o.writer->close();
        }
        
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<RExpress::Stats> stats;
            
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
