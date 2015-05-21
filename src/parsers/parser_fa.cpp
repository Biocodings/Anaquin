#include "reader.hpp"
#include "parsers/parser_fa.hpp"
#include <boost/algorithm/string.hpp>

using namespace Spike;

void ParserFA::parse(const std::string &file, std::function<void(const FALine &, const ParserProgress &)> x)
{
    Reader f(file);

    FALine l;
    std::string s;
    ParserProgress p;

    std::stringstream ss;
    #define CALL_BACK() if (p.i) { l.seq = ss.str(); x(l, p); ss.str(""); }

    while (f.nextLine(s))
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