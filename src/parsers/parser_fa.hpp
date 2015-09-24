#ifndef PARSER_FA_HPP
#define PARSER_FA_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct FALine
    {
        std::string id;
        std::string seq;
    };

    struct ParserFA
    {
        typedef std::function<void(const FALine &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif