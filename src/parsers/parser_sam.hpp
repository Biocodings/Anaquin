#ifndef PARSER_SAM_HPP
#define PARSER_SAM_HPP

#include "data/alignment.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserSAM
    {
        static void parse(const FileName &file, std::function<void (const Alignment &, const ParserProgress &)>);
    };
}

#endif