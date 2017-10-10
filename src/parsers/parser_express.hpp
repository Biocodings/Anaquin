#ifndef PARSER_EXPRESS_HPP
#define PARSER_EXPRESS_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    struct ParserExpress
    {
        static bool isExpress(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;

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
    };
}

#endif
