#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <stdexcept>

namespace Anaquin
{
    #define ASSERT(condition, message) \
    if (!(condition)) { \
        throw std::runtime_error(message); \
    }
    
    #define THROW(message) \
    throw std::runtime_error(message);
}

#endif