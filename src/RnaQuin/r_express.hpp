#ifndef R_EXPRESS_HPP
#define R_EXPRESS_HPP

#include "data/data.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    struct RExpress : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        
        enum class Format
        {
            GTF,
            Kallisto,
            Anaquin
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Format format;

            // What mixture to analyze?
            Mixture mix = Mixture::Mix_1;
        };

        struct MappingStats
        {
            // Make sure they're not duplicate
            std::set<GeneID> nGEndo;
            
            // Make sure they're not duplicate
            std::set<IsoformID> nIEndo;

            Counts nISeqs = 0;
            Counts nGSeqs = 0;
        };

        struct Stats : public RExpress::MappingStats
        {
            struct GenData
            {
                // Eg: FPKM
                double abund = NAN;
            };

            // Do we have any isoform data?
            inline bool hasIsoform() const
            {
                return isos.size() > genes.size();
            }
            
            SequinStats isos, genes;

            Limit iLimit, gLimit;
            
            // Data for the genome
            std::map<ChrID, GenData> gData;
        };

        static Stats analyze(const FileName &, const Options &o);

        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<RExpress::Stats> stats;
            
            for (const auto &file : files)
            {
                const auto x = analyze(file, o);
                
                if (x.genes.empty() && x.isos.empty())
                {
                    throw std::runtime_error("Failed to find anything on the in-silico chromosome: " + file);
                }

                stats.push_back(x);
            }

            return stats;
        }

        static Scripts generateITSV(const std::vector<RExpress::Stats> &, const Options &);
        static Scripts generateGTSV(const std::vector<RExpress::Stats> &, const Options &);

        static Scripts generateSummary(const std::vector<FileName> &,
                                       const std::vector<Stats> &,
                                       const Options &,
                                       const Units &);

        static Scripts generateRLinear(const FileName &,
                                       const std::string &,
                                       const std::vector<Stats> &,
                                       const Options &);

        static void writeITSV(const FileName &, const std::vector<RExpress::Stats> &, const Options &);
        static void writeGTSV(const FileName &, const std::vector<RExpress::Stats> &, const Options &);

        static void writeSummary(const FileName &,
                                 const std::vector<FileName> &,
                                 const std::vector<Stats> &,
                                 const Options &,
                                 const Units &);

        static void writeRLinear(const FileName &,
                                 const FileName &,
                                 const std::string &,
                                 const std::vector<Stats> &,
                                 const Options &);
        
        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
