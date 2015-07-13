#include "data/reader.hpp"
#include "parsers/parser_fa.hpp"
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

void ParserFA::parse(const Reader &r, Callback c)
{
    FALine l;
    std::string s;
    ParserProgress p;

    std::stringstream ss;
    #define CALL_BACK() if (p.i) { l.seq = ss.str(); c(l, p); ss.str(""); }

    while (r.nextLine(s))
    {
        if (s[0] != '>')
        {
            boost::trim(s);
            ss << s;
        }
        else
        {
            CALL_BACK();
            l.id = s.substr(1, s.size() - 1);
        }

        p.i++;
    }

    CALL_BACK();
}