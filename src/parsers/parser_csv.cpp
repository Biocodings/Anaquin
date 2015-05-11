#include <fstream>
#include "tokens.hpp"
#include "parser_csv.hpp"

using namespace Spike;

void ParserCSV::parse(const std::string &file, std::function<void (const std::vector<std::string> &, const ParserProgress &)> x)
{
    std::string line;
    std::ifstream i(file);
    std::vector<std::string> tokens;
    ParserProgress p;
    
    while (std::getline(i, line))
    {
        Tokens::split(line, "\t", tokens);
        x(tokens, p);
        p.i++;
    }
}