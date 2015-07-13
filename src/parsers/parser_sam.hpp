#ifndef GI_PARSER_SAM_HPP
#define GI_PARSER_SAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserSAM
    {
        static void parse(const std::string &file, std::function<void (const Alignment &, const ParserProgress &)>);
    };
}

#endif