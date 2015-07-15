#ifndef GI_PARSER_HPP
#define GI_PARSER_HPP

#include <string>
#include <stdexcept>

namespace Anaquin
{
    struct ParserProgress
    {
        long long i = 0;
    };

    struct InvalidFileError : public std::exception
    {
        InvalidFileError(const std::string &file) : file(file) {}

        const std::string file;
    };
}

#endif