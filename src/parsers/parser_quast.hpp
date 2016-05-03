#ifndef PARSER_QUAST_HPP
#define PARSER_QUAST_HPP

#include "data/types.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserQuast
    {
        struct GenomeData
        {
            // Eg: MG_29
            SequinID id;
            
            // Total length
            Base total;
            
            // Covered length
            Base covered;
        };

        /*
         *  genome_info.txt
         *
         *  reference chromosomes:
         *      CME003.v013_MG_29 (total length: 2974 bp, maximal covered length: 0 bp)
         *      CME003.v013_M3_G (total length: 1824 bp, maximal covered length: 1793 bp)
         */

        static void parseGenomeInfo(const Reader &r, std::function<void(const GenomeData &, const ParserProgress &)> f)
        {
            protectParse("genome_info.txt", [&]()
            {
                GenomeData x;
                ParserProgress p;
                
                Line line;
                std::vector<Tokens::Token> toks;
                
                while (r.nextLine(line))
                {
                    // Skip: "referebce chromosome:"
                    if (p.i++ == 0)
                    {
                        continue;
                    }
                    else if (line.find("total length") == std::string::npos)
                    {
                        continue;
                    }
                    
                    boost::trim(line);

                    /*
                     * Eg: CME003.v013_MG_29 (total length: 2974 bp, maximal covered length: 100 bp)
                     */
                    
                    Tokens::split(line, " ", toks);
                    
                    // Eg: 2974
                    x.total = stod(toks[3]);
                    
                    // Eg: 100
                    x.covered = stod(toks[8]);
                    
                    // Eg: CME003.v013_MG_29
                    Tokens::split(toks[0], "_", toks);
                    
                    // Eg: MG_29
                    x.id = toks[1] + "_" + toks[2];
                    
                    f(x, p);
                }
            });
        }
    };
}

#endif