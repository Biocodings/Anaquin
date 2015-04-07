#ifndef GI_TOKENS_HPP
#define GI_TOKENS_HPP

#include <string>
#include <boost/algorithm/string.hpp>

namespace Spike
{
    struct Tokens
    {
        template <typename T> static void split(const std::string &str, const std::string &d, T &r)
        {
            boost::split(r, str, boost::is_any_of(d));
        }
    };
}

#endif