#ifndef GI_PARSER_GTF_HPP
#define GI_PARSER_GTF_HPP

#include "feature.hpp"
#include <functional>
#include "parsers/parser.hpp"

namespace Spike
{
    struct ParserGTF
    {
        static void parse(const std::string &file, std::function<void (const Feature &)>);
    };
}

#endif