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
        
        struct Stats : public MappingStats
        {
            struct GenData
            {
                // Eg: FPKM
                double abund = NAN;
            };

            // Histogram distribution
            std::map<ChrID, Hist> isosHist, geneHist;

            LinearStats isos, genes;

            // Data for the genome
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
