#ifndef PARSER_EXPRESS_HPP
#define PARSER_EXPRESS_HPP

#include <fstream>
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

        template <typename F> static void parse(const Reader &r, F f)
        {
            ParserProgress p;
            std::string line;
            std::vector<Tokens::Token> toks;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                Data x;

                if (p.i)
                {
                    x.id    = toks[Field::Name];
                    x.cID   = toks[Field::ChrID];
                    x.abund = stod(toks[Field::Abund]);
                    
                    if (toks[Field::Start] != "-" && toks[Field::End] != "-")
                    {
                        x.l.end   = stold(toks[Field::End]);
                        x.l.start = stold(toks[Field::Start]);
                    }
                    else
                    {
                        x.l = Locus(0, 0);
                    }

                    f(x, p);
                }

                p.i++;
            }
        }
        
        static bool isIsoform(const Reader &r)
        {
            std::vector<Tokens::Token> toks;
            
            if (r.nextTokens(toks, "\t"))
            {
                return toks[1] == "IsoID";
            }
            
            return true;
        }
    };
}

#endif
