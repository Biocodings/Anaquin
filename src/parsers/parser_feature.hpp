#ifndef GI_PARSER_FEATURE_HPP
#define GI_PARSER_FEATURE_HPP

#include "parsers/parser_bed.hpp"
#include "parsers/parser_gtf.hpp"

namespace Anaquin
{
    struct ParserFeature
    {
        static void parse(const Reader &r, ParserGTF::Callback f)
        {
            ParserGTF::parse(r, f);
        }
    };
}

#endif