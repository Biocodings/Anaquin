#include <fstream>
#include <assert.h>
#include "ParserCSV.hpp"
#include <boost/algorithm/string.hpp>

bool ParserCSV::parse(const std::string &file, std::function<void (const std::vector<std::string> &)> x)
{
    std::string line;
    std::ifstream i(file);

    if (!i)
    {
        return false;
    }

    std::vector<std::string> tokens;

    while (std::getline(i, line))
    {
        boost::split(tokens, line, boost::is_any_of(","));
        x(tokens);
    }

    return true;
}