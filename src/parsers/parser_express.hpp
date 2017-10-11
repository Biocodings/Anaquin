#ifndef PARSER_EXPRESS_HPP
#define PARSER_EXPRESS_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"

namespace Anaquin
{
    struct ParserExpress
    {
        enum Field
        {
            ChrID,
            GeneID,
            IsoID,
            Abund,
        };

        struct Data
        {
            ::Anaquin::ChrID cID;
            
            // Eg: R1_1_1
            IsoformID iID;

            // Eg: R1_1
            ::Anaquin::GeneID gID;
            
            // Eg: FPKM, counts
            double abund;
        };

        static bool isExpress(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;

            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);

                if (toks.size() == 4       &&
                    toks[0] == "ChrID"     &&
                    toks[1] == "GeneID"    &&
                    toks[2] == "IsoformID" &&
                    toks[3] == "Abund")
                {
                    return true;
                }
            }
            
            return false;
        }
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            ParserProgress p;
            std::string line;
            std::vector<Token> toks;
            
            while (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                Data x;
                
                if (p.i)
                {
                    typedef ParserExpress::Field Field;

                    x.iID   = toks[Field::IsoID];
                    x.gID   = toks[Field::GeneID];
                    x.cID   = toks[Field::ChrID];
                    x.abund = ss2ld(toks[Field::Abund]);
                    f(x, p);
                }
                
                p.i++;
            }
        }
    };
}

#endif
