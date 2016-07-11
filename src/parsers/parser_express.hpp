#ifndef PARSER_EXPRESS_HPP
#define PARSER_EXPRESS_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "data/convert.hpp"

namespace Anaquin
{
    struct ParserExpress
    {
        typedef enum
        {
            ChrID,
            GeneID,
            IsoID,
            Abund,
        } Field;

        struct Data
        {
            ::Anaquin::ChrID cID;

            // Eg: R1_1 or R1_1
            GenericID id;
            
            // Eg: FPKM, counts
            double abund;
        };

        static bool isExpress(const Reader &r)
        {
            std::string line;
            std::vector<Tokens::Token> toks;

            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);

                if (toks.size() == 4 &&
                    toks[0] == "ChrID" &&
                    toks[1] == "GeneID" &&
                    toks[2] == "IsoformID" &&
                    toks[3] == "Abund")
                {
                    return true;
                }
            }
            
            return false;
        }
        
        template <typename F> static void parse(const Reader &r, bool shouldGene, F f)
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
                    if (shouldGene)
                    {
                        x.id = toks[Field::GeneID];
                    }
                    else
                    {
                        x.id = toks[Field::IsoID];
                    }
                    
                    x.cID   = toks[Field::ChrID];
                    x.abund = ns2ld(toks[Field::Abund]);
                    
                    f(x, p);
                }

                p.i++;
            }
        }
    };
}

#endif
