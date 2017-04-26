#ifndef PARSER_SALMON_HPP
#define PARSER_SALMON_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "tools/tools.hpp"
#include "data/standard.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserSalmon
    {
        enum Field
        {
            Name,
            Length,
            EffLength,
            TPM,
            NumReads,
        };

        struct Data
        {
            SequinID name;

            // Observed abundance
            Measured abund;
        };
        
        static void parse(const Reader &rr, std::function<void(const Data &, const ParserProgress &)> f)
        {
            protectParse("Salmon format", [&]()
            {
                Data d;
                ParserProgress p;
                
                Line line;
                std::vector<Token> toks;
                
                while (rr.nextLine(line))
                {
                    if (p.i++ == 0)
                    {
                        continue;
                    }
                    
                    Tokens::split(line, "\t", toks);
                    
                    d.name  = toks[Name];
                    d.abund = s2d(toks[TPM]);
                    
                    f(d, p);
                }
            });
        }
    };
}

#endif
