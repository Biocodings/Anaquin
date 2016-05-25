#ifndef PARSER_STAMP_HPP
#define PARSER_STAMP_HPP

#include "data/types.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserStamp
    {
        struct Data
        {
            GeneID id;
            
            Probability p;
            
            // Adjusted for multiple corrections
            Probability q;
            
            // Effect size
            double effect;
        };

        static void parse(const Reader &r, std::function<void(const Data &, const ParserProgress &)> f)
        {
            protectParse("STAMP TSV format", [&]()
            {
                Data d;
                ParserProgress p;
                
                std::vector<std::string> toks;
                
                while (r.nextTokens(toks, "\t"))
                {
                    if (++p.i != 1)
                    {
                        d.id = toks[0];
                        d.p  = stold(toks[1]);
                        d.q  = stold(toks[2]);
                        d.effect = stod(toks[3]);
                        
                        f(d, p);
                    }
                }
            });
        }
    };
}

#endif