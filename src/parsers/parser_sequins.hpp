#ifndef GI_PARSER_SEQUINS_HPP
#define GI_PARSER_SEQUINS_HPP

#include <string>
#include <vector>

namespace Spike
{
    typedef std::vector<std::string> SequinList;

    struct ParserSequins
    {
        static SequinList parse(const std::string &file);
    };
}

#endif