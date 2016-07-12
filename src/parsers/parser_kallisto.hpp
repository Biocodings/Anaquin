#ifndef PARSER_KALLISTO_HPP
#define PARSER_KALLISTO_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserKallisto
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
            IsoformID id;

            // Estimated abundance
            Coverage abund;
        };
        
        static bool isKallisto(const Reader &r)
        {
            std::string line;
            std::vector<Tokens::Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 5         &&
                    toks[0]  == "target_id"  &&
                    toks[1]  == "length"     &&
                    toks[2]  == "eff_length" &&
                    toks[3]  == "est_counts" &&
                    toks[4]  == "tpm")
                {
                    return true;
                }
            }

            return false;
        }

        static void parse(const Reader &r, std::function<void(const Data &, const ParserProgress &)> f)
        {
            protectParse("Kallisto TSV format", [&]()
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
                    
                    d.id    = toks[TargetID];
                    
                    // Normalized abunda (like FPKM)
                    d.abund = stod(toks[TPM]);
                    
                    f(d, p);
                }
            });
        }
    };
}

#endif