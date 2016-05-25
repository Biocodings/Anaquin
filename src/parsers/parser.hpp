#ifndef PARSER_HPP
#define PARSER_HPP

#include <string>
#include <stdexcept>

namespace Anaquin
{
    typedef std::string Line;
    
    struct ParserProgress
    {
        long long i = 0;
    };

    struct InvalidFileError : public std::exception
    {
        InvalidFileError(const std::string &file) : file(file) {}

        const std::string file;
    };
    
    template <typename F> void protectParse(const std::string &msg, F f)
    {
        try
        {
            f();
        }
        catch (...)
        {
            throw std::runtime_error("Invalid file, expected " + msg + ". Please check and try again.");
        }
    }
}

#endif