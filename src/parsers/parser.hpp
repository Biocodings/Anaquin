#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

#include <string>
#include <stdexcept>

namespace Spike
{
    struct ParserProgress
    {
        long long i = 0;
    };

    struct EmptyFileError : public std::runtime_error
    {
        EmptyFileError(const std::string &file) : std::runtime_error(file) {}
    };

    struct InvalidFileError : public std::runtime_error
    {
        InvalidFileError(const std::string &file) : std::runtime_error(file) {}
    };

    enum DataMode
    {
        File,
        String,
    };
}

#endif