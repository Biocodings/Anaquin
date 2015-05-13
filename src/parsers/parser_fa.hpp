#ifndef GI_PARSER_FA_HPP
#define GI_PARSER_FA_HPP

#include "analyzer.hpp"

namespace Spike
{
    struct FALine
    {
        std::string id;
        std::string seq;
    };

    struct ParserFA
    {
        static void parse(const std::string &file, std::function<void(const FALine &, const ParserProgress &)>);
    };
}

#endif