#ifndef PARSER_KALLISTO_HPP
#define PARSER_KALLISTO_HPP

#include "data/types.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParseKallisto
    {
        enum Field
        {
            TargetID,
            Length,
            EffLength,
            EstCounts,
            TPM,
        };

        struct Data
        {
            // Eg: R1_101_1
            IsoformID iID;

            // Estimated abundance
            Coverage cov;
        };
        
        static void parse(const Reader &r, std::function<void(const Data &, const ParserProgress &)> f)
        {
            Data d;
            ParserProgress p;
            
            Line line;
            std::vector<Tokens::Token> toks;

            while (r.nextLine(line))
            {
                if (p.i++ == 0)
                {
                    continue;
                }
                
                Tokens::split(line, "\t", toks);
                
                d.iID = toks[TargetID];
                d.cov = stod(toks[EstCounts]);

                f(d, p);
            }
        }
    };
}

#endif