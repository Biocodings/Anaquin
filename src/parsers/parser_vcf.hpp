#ifndef GI_PARSER_VCF_HPP
#define GI_PARSER_VCF_HPP

#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Spike
{
    typedef Variation VCFVariant;

    struct ParserVCF
    {
        typedef std::function<void (const VCFVariant &, const ParserProgress &)> Callback;
        static void parse(const Reader &r, Callback);
    };
}

#endif