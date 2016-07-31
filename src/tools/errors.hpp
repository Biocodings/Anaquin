#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <stdexcept>

namespace Anaquin
{
    #define A_ASSERT(condition, message) \
    if (!(condition)) { \
        throw std::runtime_error(message); \
    }
    
    #define A_THROW(message) \
    throw std::runtime_error(message);
}

#endif