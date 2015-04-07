#ifndef GI_PARSER_SAM_HPP
#define GI_PARSER_SAM_HPP

#include <functional>
#include "alignment.hpp"

namespace Spike
{
    struct ParserSAM
    {
        static void parse(const std::string &file, std::function<void (const Alignment &)>);
    };
}

#endif