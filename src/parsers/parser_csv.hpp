#ifndef GI_PARSER_CSV_HPP
#define GI_PARSER_CSV_HPP

#include <vector>
#include <functional>
#include "parsers/parser.hpp"

namespace Spike
{
    typedef std::vector<std::string> Fields;

    struct ParserCSV
    {
        typedef std::function<void (const Fields &, const ParserProgress &)> Callback;
        static void parse(const std::string &file, Callback, DataMode mode, const std::string &delim = "\t");
    };    
}

#endif