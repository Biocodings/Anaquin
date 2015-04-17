#ifndef GI_INVALID_FILE_HPP
#define GI_INVALID_FILE_HPP

#include <string>
#include <stdexcept>

namespace Spike
{
    class InvalidFile : public std::exception
    {
        public:
            InvalidFile(const std::string &file) : std::exception() {}
    };
}

#endif