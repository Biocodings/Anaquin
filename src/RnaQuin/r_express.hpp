/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef R_EXPRESS_HPP
#define R_EXPRESS_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser_cufflink.hpp"

namespace Anaquin
{
    struct RExpress : public Analyzer
    {
        typedef ParserCufflink::Data TestData;
        
        enum class Metrics
        {
            Gene,
            Isoform,
        };
        
        enum class Format
        {
            GTF,
            Kallisto,
            Text
        };
        
        struct Options : public AnalyzerOptions
        {
            Options() {}

            Format format;

            // What mixture to analyze?
            Mixture mix = Mixture::Mix_1;

            // Gene or isoform?
            Metrics metrs;
        };
        
        struct Stats : public MappingStats, public LimitStats
        {
            struct GenData
            {
                // Eg: FPKM
                double abund = NAN;
            };

            SequinStats isos, genes;

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
                
                if (x.genes.empty() && x.isos.empty() && files.size() == 1)
                {
                    throw std::runtime_error("Failed to find anything on the in-silico chromosome: " + file);
                }

                stats.push_back(x);
            }

            return stats;
        }

        static Scripts generateCSV(const std::vector<RExpress::Stats> &,
                                   const Options &);

        static Scripts generateSummary(const std::vector<FileName> &,
                                       const std::vector<Stats> &,
                                       const Options &,
                                       const Units &);

        static Scripts generateRLinear(const FileName &,
                                       const std::vector<Stats> &,
                                       const Options &);

        static void writeCSV(const FileName &, const std::vector<RExpress::Stats> &, const Options &);

        static void writeSummary(const FileName &,
                                 const std::vector<FileName> &,
                                 const std::vector<Stats> &,
                                 const Options &,
                                 const Units &);

        static void writeRLinear(const FileName &,
                                 const FileName &,
                                 const std::vector<Stats> &,
                                 const Options &);
        
        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
