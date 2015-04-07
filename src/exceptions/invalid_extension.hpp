#ifndef GI_INVALID_EXTENSION_HPP
#define GI_INVALID_EXTENSION_HPP

#include <string>
#include <stdexcept>

namespace Spike
{
    class InvalidExtension : public std::exception
    {
        public:
            InvalidExtension(const std::string &file, const std::string &expected, const std::string &actual)
                : std::exception(), file(file), expected(expected), actual(actual) {}

        private:
            const std::string file, expected, actual;
    };
}


#endif