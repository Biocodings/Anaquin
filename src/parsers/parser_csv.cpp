#include "parsers/parser_csv.hpp"

using namespace Anaquin;

void ParserCSV::parse(const Reader &r, Callback c, const std::string &delim)
{
    ParserProgress p;
    std::vector<std::string> tokens;

    while (r.nextTokens(tokens, delim))
    {
        c(tokens, p);
        p.i++;
    }
}