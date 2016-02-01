#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include "data/variant.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef Variant VCFVariant;
        typedef std::function<void (const VCFVariant &, const ParserProgress &)> Callback;

        static void parse(const Reader &, Callback);
    };
}

#endif