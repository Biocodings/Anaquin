#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <stdexcept>

namespace SS
{
    #define ASSERT(cond, message) \
        if (!(cond)) { \
            throw std::runtime_error(message); \
        }
}

#endif