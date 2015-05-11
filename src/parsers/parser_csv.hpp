#ifndef GI_PARSER_CSV_HPP
#define GI_PARSER_CSV_HPP

#include <functional>
#include "analyzer.hpp"

namespace Spike
{
    typedef std::vector<std::string> Fields;

    struct ParserCSV
    {
        static void parse(const std::string &file, std::function<void (const Fields &, const ParserProgress &)>);
    };    
}

#endif