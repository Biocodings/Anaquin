#include "reader.hpp"
#include "parser_csv.hpp"

using namespace Spike;

void ParserCSV::parse(const std::string &file, Callback c, ParserMode mode, const std::string &delim)
{
    ParserProgress p;
    Reader r(file, mode);
    std::vector<std::string> tokens;
    
    while (r.nextTokens(tokens, delim))
    {
        c(tokens, p);
        p.i++;
    }
}