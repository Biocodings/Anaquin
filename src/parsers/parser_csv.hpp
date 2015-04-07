#ifndef GI_PARSER_CSV_HPP
#define GI_PARSER_CSV_HPP

#include <vector>
#include <functional>

namespace Spike
{
    typedef std::vector<std::string> Fields;
    
    struct ParserCSV
    {
        static bool parse(const std::string &file, std::function<void (const Fields &)>);
    };    
}

#endif