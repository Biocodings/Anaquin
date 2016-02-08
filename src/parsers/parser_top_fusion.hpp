#ifndef PARSER_TOP_FUSION_HPP
#define PARSER_TOP_FUSION_HPP

#include <functional>
#include "data/biology.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserTopFusion
    {
        typedef CalledFusion Data;
        typedef std::function<void(const Data &, const ParserProgress &)> Callback;

        static void parse(const Reader &, Callback);
    };
}

#endif