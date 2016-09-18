#ifndef PARSER_SALMON_HPP
#define PARSER_SALMON_HPP

#include "data/data.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/convert.hpp"
#include "data/standard.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserSalmon
    {
        enum Field
        {
            Name,
            Length,
            EffLength,
            TPM,
            NumReads,
        };

        struct Data
        {
            ChrID cID;
            
            SequinID id;

            // Estimated abundance
            Measured abund;
        };
        
        static bool isSalmon(const Reader &r)
        {
            std::string line;
            std::vector<Token> toks;
            
            // Read the header
            if (r.nextLine(line))
            {
                Tokens::split(line, "\t", toks);
                
                if (toks.size() == 5              &&
                    toks[0]  == "Name"            &&
                    toks[1]  == "Length"          &&
                    toks[2]  == "EffectiveLength" &&
                    toks[3]  == "TPM"             &&
                    toks[4]  == "NumReads")
                {
                    return true;
                }
            }

            return false;
        }

        static void parse(const Reader &rr, std::function<void(const Data &, const ParserProgress &)> f)
        {
            protectParse("Salmon format", [&]()
            {
                //const auto &r = Standard::instance().v_ref;

                Data d;
                ParserProgress p;
                
                Line line;
                std::vector<Token> toks;
                
                while (rr.nextLine(line))
                {
                    if (p.i++ == 0)
                    {
                        continue;
                    }
                    
                    Tokens::split(line, "\t", toks);
                    
                    d.id = toks[Name];
                    
                    //if (r.findTrans(ChrIS, d.id))
                    {
                        d.cID = ChrIS;
                    }
                    //else
                    //{
                        // We don't know exactly where it is...
                    //    d.cID = Geno;
                    //}

                    d.abund = s2d(toks[TPM]);
                    
                    f(d, p);
                }
            });
        }
    };
}

#endif