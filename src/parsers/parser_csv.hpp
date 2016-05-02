#ifndef PARSER_CSV_HPP
#define PARSER_CSV_HPP

#include <vector>
#include <functional>
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserCSV
    {
        typedef std::vector<std::string> Data;

        // Callback for parsing a CSV file
        typedef std::function<void (const ParserCSV::Data &, const ParserProgress &)> Callback;

        static void parse(const Reader &r, Callback c, const std::string &delim = ",")
        {
            protectParse("CSV", [&]()
            {
                ParserProgress p;
                std::vector<std::string> tokens;
                
                while (r.nextTokens(tokens, delim))
                {
                    c(tokens, p);
                    p.i++;
                }
            });
        }
    };
}

#endif