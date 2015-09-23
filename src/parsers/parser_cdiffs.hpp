#ifndef GI_PARSER_CDIFFS_HPP
#define GI_PARSER_CDIFFS_HPP

#include <ss/p.hpp>
#include "stats/analyzer.hpp"
#include "parsers/parser_cuffs.hpp"

namespace Anaquin
{
    struct TrackingDiffs
    {
        TrackID chromID;
        TrackID testID;
        GeneID  geneID;

        FPKM fpkm_1, fpkm_2;

        SS::LogFold logFold;
        SS::TestStats stats;

        // The p-value and q-value under the null-hypothesis
        SS::P p, q;

        TrackingStatus status;
    };

    struct ParserCDiffs
    {
        static void parse(const FileName &file, std::function<void (const TrackingDiffs &, const ParserProgress &)>);
    };    
}

#endif
