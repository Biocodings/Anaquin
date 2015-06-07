#ifndef GI_PARSER_TMAP_HPP
#define GI_PARSER_TMAP_HPP

#include <functional>
#include "stats/analyzer.hpp"

namespace Spike
{
    typedef std::string RefGeneID;
    typedef std::string RefID;

    struct TMap
    {
        RefGeneID refGenID;
        RefID refID;

        FPKM fpkm;
        FPKM lFPKM;
        FPKM uFPKM;
    };

    struct ParserTMap
    {
        static void parse(const std::string &file, std::function<void (const TMap &, const ParserProgress &)>);
    };
}

#endif