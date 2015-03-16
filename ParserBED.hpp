#ifndef AS_PARSER_BED_HPP
#define AS_PARSER_BED_HPP

#include "Feature.hpp"
#include <functional>

struct ParserBED
{
    static bool parse(const std::string &file, std::function<void(const Feature &)>);
};

#endif