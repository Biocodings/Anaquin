#ifndef PARSER_EXPRESS_HPP
#define PARSER_EXPRESS_HPP

#include "data/types.hpp"
#include "data/tokens.hpp"

namespace Anaquin
{
    struct ParserExpress
    {
        typedef enum
        {
            ChrID,
            Name,
            Start,
            End,
            Abund,
        } Field;

        struct Data
        {
            ::Anaquin::ChrID cID;

            // Eg: genes or isoforms
            std::string id;
            
            // Position of the gene/isoform
            Locus l;
            
            // Eg: FPKM, counts
            double abund;
        };

        template <typename F> static void parse(const FileName &file, F f)
        {
            Reader r(file);
            ParserProgress p;
            std::vector<Tokens::Token> toks;
            std::string line;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                Data x;

                if (p.i)
                {
                    x.cID     = toks[Field::ChrID];
                    x.id      = toks[Field::Name];
                    x.l.start = stold(toks[Field::Start]);
                    x.l.end   = stold(toks[Field::End]);
                    x.abund   = stod(toks[Field::Abund]);

                    f(x, p);
                }

                p.i++;
            }
        }
    };
}

#endif
