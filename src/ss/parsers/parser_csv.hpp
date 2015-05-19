#ifndef SS_PARSER_CSV_HPP
#define SS_PARSER_CSV_HPP

#include <string>
#include <vector>
#include <fstream>
#include <boost/algorithm/string.hpp>

namespace SS
{
    struct ParserCSV
    {
        typedef std::vector<std::string> CSVFields;

        template<typename F> static void parse(const std::string &file, F f)
        {
            std::string l;
            std::ifstream i(file);
            std::vector<std::string> toks;

            while (std::getline(i, l))
            {
                toks.clear();
                boost::split(toks, l, boost::is_any_of(","));
                f(toks);
            }
        }
    };    
}

#endif