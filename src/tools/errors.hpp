#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <assert.h>
#include <stdexcept>

namespace Anaquin
{
    #define A_ASSERT(cond) \
    assert(cond);

    #define A_CHECK(cond, message) \
    if (!(cond)) { throw std::runtime_error(message); }

    #define A_THROW(message) \
    throw std::runtime_error(message);
}

#endif