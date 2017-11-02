#ifndef PARSER_FA_HPP
#define PARSER_FA_HPP

#include "data/reader.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct ParserFA
    {
        struct Data
        {
            ChrID id;
            Sequence seq;
        };

        typedef std::function<void(const Data &, const ParserProgress &)> Callback;

        static void parse(const Reader &r, Callback f, const ChrID &chrID = "")
        {
            Data l;
            std::string s;
            ParserProgress p;
            
            std::stringstream ss;
            #define CALL_BACK() if (p.i) { l.seq = ss.str(); f(l, p); ss.str(""); }

            while (r.nextLine(s))
            {
                if (s[0] != '>')
                {
                    if (chrID.empty() || l.id == chrID)
                    {
                        boost::trim(s);
                        ss << s;
                    }
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
    };
}

#endif
