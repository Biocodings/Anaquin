#include <assert.h>
#include "file.hpp"
#include "parsers/parser_fa.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

void ParserFA::parse(const std::string &file, std::function<void(const FALine &, const ParserProgress &)> x)
{
    File f(file);
    FALine l;
    ParserProgress p;
    std::vector<std::string> tokens;
    
    while (f.nextTokens(tokens, "\t"))
    {
        if (p.i % 2)
        {
            x(l, p);
        }
        else
        {
            l.id = tokens[0];
        }

        p.i++;
    }
}