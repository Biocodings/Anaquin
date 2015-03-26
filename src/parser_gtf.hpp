#ifndef AS_PARSER_GTF_HPP
#define AS_PARSER_GTF_HPP

#include "Feature.hpp"
#include <functional>

struct ParserProgress
{
    Lines i = 0;
    
    // Whether parsing should be stopped
    bool terminate = false;
};

struct ParserGTF
{
    static bool parse(const std::string &file, std::function<void (const Feature &, ParserProgress &p)>);
};

#endif