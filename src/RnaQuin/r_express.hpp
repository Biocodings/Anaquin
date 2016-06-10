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
extern Anaquin::Scripts PlotScatter();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMultiple();

// Defined in resources.cpp
extern Anaquin::Scripts PlotTMinor();

extern Anaquin::FileName MixRef();

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
            
            // What software generated the files?
            Software soft;
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
            
            if (stats.size() == 1)
            {
                o.writer->write(RWriter::createScript(csv, PlotScatter()));
            }
            else
            {
                o.writer->write(RWriter::createScript(csv, PlotTMultiple()));
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
            o.info("Generating " + summary);
            o.writer->open(summary);
            
            const auto format = "-------RnaExpression Output\n\n"
                                "Summary for input: %1%\n\n"
                                "*Arithmetic average and standard deviation are shown\n\n"
                                "-------User Transcript Annotations\n\n"
                                "Annotation file: %2%\n"
                                "Synthetic: %3%\n"
                                "Genome:    %4%\n\n"
                                "Mixture file: %5%\n\n"
                                "-------Genes Expressed\n\n"
                                "       Synthetic: %6%\n"
                                "       Detection Sensitivity: %7% (attomol/ul) (%8%)\n\n"
                                "       Genome: %9%\n\n"
                                "-------Limit of Quantification (LOQ)\n"
                                "       *Estimated by piecewise segmented regression\n\n"
                                "Break: %10% (%11%)\n\n"
                                "*Below LOQ\n"
                                "Intercept:   %12%\n"
                                "Slope:       %13%\n"
                                "Correlation: %14%\n"
                                "R2:          %15%\n\n"
                                "*Above LOQ\n"
                                "Intercept:   %16%\n"
                                "Slope:       %17%\n"
                                "Correlation: %18%\n"
                                "R2:          %19%\n\n"
                                "-------Linear regression (log2 scale)\n\n"
                                "Correlation: %20%\n"
                                "Slope:       %21%\n"
                                "R2:          %22%\n"
                                "F-statistic: %23%\n"
                                "P-value:     %24%\n"
                                "SSM:         %25%, DF: %26%\n"
                                "SSE:         %27%, DF: %28%\n"
                                "SST:         %29%, DF: %30%\n";
        
            o.writer->write((boost::format(format) % "????"
                                                   % "????"
                                                   % "????"   // 3
                                                   % "????"   // 4
                                                   % MixRef() // 5
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????" // 10
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????" // 15
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????" // 20
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????" // 25
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????"
                                                   % "????" // 30
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
