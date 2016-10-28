#ifndef PARSER_HPP
#define PARSER_HPP

#include "tools/errors.hpp"

namespace Anaquin
{
    typedef std::string Line;
    
    struct ParserProgress
    {
        // Whether parsing should be stopped immediately
        bool stopped = false;
        
        long long i = 0;
    };

    template <typename F> void protectParse(const std::string &msg, F f)
    {
        try
        {
            f();
        }
        catch (...)
        {
            throw InvalidFormatException("Invalid file, expected " + msg);
        }
    }
}

#endif
