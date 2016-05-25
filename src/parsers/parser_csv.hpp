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

        template <typename F> static void parse(const Reader &r, F f, const std::string &delim = ",")
        {
            protectParse("CSV format", [&]()
            {
                ParserProgress p;
                std::vector<std::string> tokens;
                
                while (r.nextTokens(tokens, delim))
                {
                    f(tokens, p);
                    p.i++;
                }
            });
        }
    };
}

#endif