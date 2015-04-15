#ifndef GI_PARSER_CDIFFS_HPP
#define GI_PARSER_CDIFFS_HPP

#include "types.hpp"
#include <functional>
#include <ss/prob.hpp>
#include "parsers/parser_cuffs.hpp"

namespace Spike
{
    struct TrackingDiffs
    {
        TrackID testID;
        GeneID  geneID;

        FPKM fpkm_1, fpkm_2;

        SS::LogFold logFold;
        SS::TestStats stats;

        // The p-value and q-value under the null-hypothesis
        SS::Probability p, q;

        TrackingStatus status;
    };

    struct ParserCDiffs
    {
        static void parse(const std::string &file, std::function<void (const TrackingDiffs &)>);
    };    
}

#endif
