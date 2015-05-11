#include <assert.h>
#include "d_align.hpp"
#include "parsers/parser_sam.hpp"
#include <iostream>
using namespace Spike;

DAlignStats DAlign::analyze(const std::string &file, const Options &options)
{
    DAlignStats stats;
    const auto &s = Standard::instance();

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        if (classify(stats.m, align, [&](const Alignment &)
                     {
                         const auto matched = find(s.d_exons, align, MatchRule::Contains);
                        
                         if (!matched)
                         {
                             return Negative;
                         }

                         return Positive;
                     }))
        {
            // Empty Implementation
        }
    });

    //AnalyzeReporter::report("dalign_base.stats", stats.m, stats.s, stats.c, options.writer);

	return DAlignStats();
}