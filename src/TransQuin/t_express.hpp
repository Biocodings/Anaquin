#ifndef T_REPLICATE_HPP
#define T_REPLICATE_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

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
