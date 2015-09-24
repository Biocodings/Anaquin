#ifndef PARSER_TMAP_HPP
#define PARSER_TMAP_HPP

#include <functional>
#include "stats/analyzer.hpp"

namespace Anaquin
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
        static void parse(const std::string &, std::function<void (const TMap &, const ParserProgress &)>);
    };
}

#endif