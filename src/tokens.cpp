#include "tokens.hpp"
#include <boost/algorithm/string.hpp>

void Tokens::split(const std::string &s, const std::string &d, std::vector<std::string> &ts)
{
    boost::split(ts, s, boost::is_any_of(d));
}