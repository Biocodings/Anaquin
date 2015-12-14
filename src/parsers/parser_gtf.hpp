#ifndef PARSER_GTF_HPP
#define PARSER_GTF_HPP

#include "data/feature.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserGTF
    {
        typedef std::function<void (const Feature &, const std::string &, const ParserProgress &)> Callback;
        static void parse(const Reader &, Callback);
    };
}

#endif