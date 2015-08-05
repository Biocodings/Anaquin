#ifndef GI_PARSER_CSV_HPP
#define GI_PARSER_CSV_HPP

#include <vector>
#include <functional>
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserCSV
    {
        typedef std::vector<std::string> Fields;

        // Callback for parsing a CSV file
        typedef std::function<void (const ParserCSV::Fields &, const ParserProgress &)> Callback;

        static void parse(const Reader &, Callback, const std::string &delim = ",");
    };
}

#endif