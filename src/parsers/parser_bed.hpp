#ifndef PARSER_BED_HPP
#define PARSER_BED_HPP

#include <functional>
#include "data/types.hpp"
#include "data/locus.hpp"
#include "data/reader.hpp"
#include "data/biology.hpp"
#include "parsers/parser.hpp"
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct ParserBed
    {
        struct Data
        {
            operator const std::string &() const { return name; }
            
            ChrID cID;
            
            // Forward or reverse strand?
            Strand strand;

            Locus l;

            // Eg: chr1_10482481_10483779
            std::string name;
        };

        static void parse(const Reader &r, std::function<void(const Data &, const ParserProgress &)> f)
        {
            protectParse("BED", [&]()
            {
                Data d;
                ParserProgress p;
                
                std::vector<std::string> sizes, starts, tokens;
                
                while (r.nextTokens(tokens, "\t"))
                {
                    // Empty line?
                    if (tokens.size() == 1)
                    {
                        return;
                    }
                    
                    // Name of the chromosome
                    d.cID = tokens[0];
                    
                    // Position of the feature in standard chromosomal coordinates
                    d.l = Locus(stod(tokens[1]) + 1, stod(tokens[2]));
                    
                    if (tokens.size() >= 6)
                    {
                        // Defines the strand, either '+' or '-'
                        d.strand = tokens[5] == "+" ? Forward : Backward;
                    }
                    
                    if (tokens.size() >= 4)
                    {
                        // Name of the BED line (eg: gene)
                        d.name = tokens[3];
                    }
                    
                    f(d, p);
                    p.i++;
                }
            });
        }
    };
}

#endif