#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include "tools/tools.hpp"
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/biology.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef Variant Data;

        enum Field
        {
            Chrom,
            Pos,
            ID,
            Ref,
            Alt,
            Qual,
            Filter,
            Info,
            Format,
            FormatData
        };

        template <typename F> static void parse(const Reader &r, F f)
        {
            std::string line;
            
            ParserProgress p;
            
            std::vector<std::string> tmp;
            std::vector<std::string> fields;
            std::vector<std::string> formats;
            
            while (r.nextLine(line))
            {
                Data x;
                p.i++;
                
                if (p.stopped)
                {
                    break;
                }                
                else if (line[0] == '#')
                {
                    continue;
                }
                
                Tokens::split(line, "\t", fields);
                
                x.cID = standChr(fields[Field::Chrom]);
                
                // Eg: D_1_3_R
                x.name = fields[Field::ID];
                
                // VCF has 1-based position
                x.l.start = x.l.end = stod(fields[Field::Pos]);
                
                // Reference allele
                x.ref = fields[Field::Ref];
                
                /*
                 * Additional information
                 *
                 *    AF: allele frequency for each ALT allele in the same order as listed
                 */
                
                if (fields[Field::Info] != ".")
                {
                    static std::vector<std::string> infos;
                    Tokens::split(fields[Field::Info], ";", infos);
                    
                    for (const auto &info : infos)
                    {
                        Tokens::split(info, "=", tmp);
                        
                        // Measured allele frequency
                        if (tmp[0] == "AF") { x.allF = stof(tmp[1]); }

                        else if (tmp.size() == 2)
                        {
                            x.opts[tmp[0]] = tmp[1];
                        }
                    }
                }

                std::vector<Sequence> alts;
                Tokens::split(fields[Field::Alt], ",", alts);
                
                // Ignore anything that is not really a variant
                if ((x.alt = alts[0]) == ".")
                {
                    continue;
                }
                
                x.qual = fields[Field::Qual] != "." ? s2d(fields[Field::Qual]) : NAN;

                if (fields.size() > Field::Format)
                {
                    Tokens::split(fields[Field::Format], ":", formats);
                    
                    // Eg: 1/2:0,11,5:16:99:694,166,119,378,0,331
                    Tokens::split(fields[Field::FormatData], ":", tmp);
                    
                    // Check all the format data...
                    for (auto j = 0; j < tmp.size(); j++)
                    {
                        if (formats[j] == "AD")
                        {
                            std::vector<std::string> toks;
                            Tokens::split(tmp[j], ",", toks);
                            
                            if (toks.size() == 1)
                            {
                                x.readR = s2d(toks[0]);
                                x.readV = s2d(toks[0]);
                            }
                            else
                            {
                                x.readR = s2d(toks[0]);
                                x.readV = s2d(toks[1]);
                            }
                        }
                        else if (formats[j] == "DP")
                        {
                            x.depth = s2d(tmp[j]);
                        }
                    }
                }
                
                f(x, p);
            }
        }
    };
}

#endif
