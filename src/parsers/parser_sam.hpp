#ifndef GI_PARSER_SAM_HPP
#define GI_PARSER_SAM_HPP

#include "analyzer.hpp"
#include "alignment.hpp"

namespace Spike
{
    struct ParserSAM
    {
        static void parse(const std::string &file, std::function<void (const Alignment &, const ParserProgress &)>);
    };
}

#endif