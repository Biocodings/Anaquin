#ifndef GI_TOKENS_HPP
#define GI_TOKENS_HPP

#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct Tokens
    {
        typedef std::string Token;
        
        template <typename T> static void split(const std::string &str, const std::string &d, T &r)
        {
            r.clear();
            boost::split(r, str, boost::is_any_of(d));
        }

        static Token first(const std::string &str, const std::string &d)
        {
            static std::vector<std::string> toks;

            toks.clear();
            Tokens::split(str, d, toks);

            return toks[0];
        }
    };
}

#endif