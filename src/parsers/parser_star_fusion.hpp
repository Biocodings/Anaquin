#ifndef PARSER_STAR_FUSION_HPP
#define PARSER_STAR_FUSION_HPP

#include <functional>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserStarFusion
    {
        typedef CalledFusion Data;
        typedef std::function<void (const Data &, const ParserProgress &)> Functor;

        // Parse an output file from FusionStar
        static void parse(const Reader &, Functor);
    };
}

#endif