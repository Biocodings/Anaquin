/*
 * Copyright (C) 2016 - Garvan Institute of Medical Research
 *
 *  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
 */

#ifndef T_DIFFS_HPP
#define T_DIFFS_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct TDiffs : public Analyzer
    {
        enum class Software
        {
            Cuffdiffs,
        };
        
        enum class Metrics
        {
            Gene,
            Isoform
        };

        struct Options : public DoubleMixtureOptions
        {
            Options() {}

            Software soft = Software::Cuffdiffs;
            Metrics metrs = Metrics::Gene;
        };

        struct Stats : public MappingStats
        {
            typedef LinearStats Data;
            
            std::map<ChromoID, Data> data;

            Limit limit;
        };

        static Stats analyze(const FileName &, const Options &o);
        static Stats analyze(const std::vector<DiffTest> &, const Options &o);

        // Analyze multiple replicates
        static std::vector<Stats> analyze(const std::vector<FileName> &files, const Options &o)
        {
            std::vector<TDiffs::Stats> stats;
            
            for (const auto &file : files)
            {
                stats.push_back(analyze(file, o));
            }
            
            return stats;            
        }

        static void report(const FileName &, const Options &o = Options());
        static void report(const std::vector<FileName> &, const Options &o = Options());
    };
}

#endif
