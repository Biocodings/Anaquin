#ifndef GI_TOKENS_HPP
#define GI_TOKENS_HPP

#include <string>
#include <vector>

struct Tokens
{
    static void split(const std::string &, const std::string &, std::vector<std::string> &);
};

#endif