#ifndef PARSER_SEQUINS_HPP
#define PARSER_SEQUINS_HPP

#include <string>
#include <vector>

namespace Anaquin
{
    typedef std::vector<std::string> SequinList;

    struct ParserSequins
    {
        static SequinList parse(const FileName &file);
    };
}

#endif