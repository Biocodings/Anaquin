#ifndef GI_PARSER_GTF_HPP
#define GI_PARSER_GTF_HPP

#include "feature.hpp"
#include <functional>

namespace Spike
{
    struct ParserProgress
    {
        Lines i = 0;
    };
    
    struct ParserGTF
    {
        static void parse(const std::string &file, std::function<void (const Feature &, ParserProgress &p)>);
    };
}

#endif