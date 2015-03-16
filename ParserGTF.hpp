#ifndef AS_PARSER_GTF_HPP
#define AS_PARSER_GTF_HPP

#include "Feature.hpp"
#include <functional>

struct ParserGTF
{
    static bool parse(const std::string &file, std::function<void(const Feature &)>);
};

#endif