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
            Isoform
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

            LinearStats isos, genes;

            // Data for the genome
            std::map<GenoID, GenData> gData;
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

        static std::vector<Stats> analyze(const std::vector<TestData> &, const Options &);
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
