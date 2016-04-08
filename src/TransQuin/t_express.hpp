#ifndef T_REPLICATE_HPP
#define T_REPLICATE_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    // Defined in resources.cpp
    extern Scripts PlotTAbundAbund();
    
    // Defined in resources.cpp
    extern Scripts PlotRAbundAbund();

    // Defined in resources.cpp
    extern Scripts PlotMajor();

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
            return StatsWriter::inflectSummary(o.rChrT,
                                               o.rEndo,
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
            
            return StatsWriter::inflectSummary(o.rChrT, o.rEndo, files, sHists, mStats, lStats, units);
        }
        
        template <typename Stats> static Scripts multipleCSV(const std::vector<Stats> &stats)
        {
            std::set<SequinID> seqs;
            
            // This is the data structure that will be convenient
            std::map<unsigned, std::map<SequinID, Concent>> data;
            
            // Expected concentration
            std::map<SequinID, Concent> expect;
            
            std::stringstream ss;
            ss << "Sequin\tEAbund";
            
            for (auto i = 0; i < stats.size(); i++)
            {
                ss << ((boost::format("\tA%1%") % (i+1)).str());
                
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

        template <typename Options> static void generateRMajor(const FileName &output,
                                                               const FileName &csv,
                                                               const std::vector<Stats> &stats,
                                                               const Options &o)
        {
            if (stats.size() >= 2 && o.metrs == TExpress::Metrics::Isoform)
            {
                o.writer->open(output);
                o.writer->write(RWriter::createScript(csv, PlotMajor()));
                o.writer->close();
            }
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
                o.writer->write(RWriter::createScript(csv, PlotTAbundAbund()));
            }
            else
            {
                o.writer->write(RWriter::createScript(csv, PlotRAbundAbund()));
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
                o.writer->write(StatsWriter::writeCSV(stats[0], "EAbund", "MAbund"));
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
            
            if (files.size() == 1)
            {
                o.writer->write(singleSummary(stats[0], files[0], units, o));
            }
            else
            {
                o.writer->write(multipleSummary(files, stats, units, o));
            }
            
            o.writer->close();
        }
        
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            if (files.size() < 2)
            {
                throw std::runtime_error("Two or more replicates required");
            }
            
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
