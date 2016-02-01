#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    typedef Variation VCFVariant;

    struct ParserVCF
    {
        typedef std::function<void (const VCFVariant &, const ParserProgress &)> Callback;

        static void parse(const Reader &, Callback);
    };
}

#endif