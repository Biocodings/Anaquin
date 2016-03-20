#ifndef PARSER_CNVNATOR_HPP
#define PARSER_CNVNATOR_HPP

#include "data/tokens.hpp"
#include "stats/analyzer.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserCNV
    {
        enum VCFField
        {
            Type,
            Location,
            Size,
            NormalizedDepth
        };
        
        struct Data
        {
            ChrID cID;
            
            // Location of the segment
            Locus l;
            
            // Size of the segment
            Base size;
            
            double fold;
        };
        
        template <typename F> static void parse(const Reader &r, F f)
        {
            /*
             * Eg: duplication	chrR:24900001-24900600	600	19.6757	0.00752628	0	1	1	1
             */
            
            std::string line;
            std::vector<std::string> toks;
            
            Data c;
            ParserProgress p;
            
            while (r.nextLine(line))
            {
                p.i++;
                Tokens::split(line, "\t", toks);
                
                if (toks.size() >= 9 && toks[0] == "duplication")
                {
                    std::vector<std::string> tmp;
                    
                    // Eg: chrR:24900001-24900600
                    Tokens::split(toks[Location], ":", tmp);
                    
                    if (tmp.size() != 2)
                    {
                        continue;
                    }
                    
                    // Eg: chrR
                    c.cID = tmp[0];
                    
                    // Eg: 24900001-24900600
                    Tokens::split(tmp[1], "-", tmp);
                    
                    if (tmp.size() != 2)
                    {
                        continue;
                    }

                    c.l = Locus(stod(tmp[0]), stod(tmp[1]));
                    
                    // Eg: 600
                    c.size = stod(toks[Size]);
                    
                    // Eg: 19.6757
                    c.fold = stod(toks[NormalizedDepth]);
                    
                    f(c, p);
                }
            }
        }
    };
}

#endif