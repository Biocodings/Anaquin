#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef CalledVariant Data;

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
            
            Data d;
            ParserProgress p;
            
            std::vector<std::string> t;
            std::vector<std::string> infos;
            std::vector<std::string> fields;
            std::vector<std::string> formats;
            
            while (r.nextLine(line))
            {
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
                
                // Eg: chrT
                d.cID = fields[Field::Chrom];
                
                // Eg: D_1_3_R
                d.id = fields[Field::ID];
                
                // VCF has 1-based position
                d.l.start = d.l.end = stod(fields[Field::Pos]);
                
                // Reference allele
                d.ref = fields[Field::Ref];
                
                /*
                 * Additional information
                 *
                 *    AF: allele frequency for each ALT allele in the same order as listed
                 */
                
                if (fields[Field::Info] != ".")
                {
                    Tokens::split(fields[Field::Info], ";", infos);
                    
                    for (const auto &info : infos)
                    {
                        /*
                         * Eg: AA=.;DP=124
                         *     AA=g;DP=132;HM2
                         */
                        
                        Tokens::split(info, "=", t);
                        
                        // Measured allele frequency
                        if (t[0] == "AF") { d.allF = stof(t[1]); }
                    }
                }
                
                Tokens::split(fields[Field::Format], ":", formats);
                assert(std::find(formats.begin(), formats.end(), "GT") != formats.end());
                
                /*
                 * Anaquin doesn't really support multi-alleles (because VarQuin design doesn't have it).
                 */
                
                std::vector<Sequence> alts;
                Tokens::split(fields[Field::Alt], ",", alts);
                
                // Simply ignore any other allele
                d.alt = alts[0];

                // Eg: 1/2:0,11,5:16:99:694,166,119,378,0,331
                Tokens::split(fields[Field::FormatData], ":", t);

                // Check all the format data...
                for (auto j = 0; j < t.size(); j++)
                {
                    if (formats[j] == "AD")
                    {
                        std::vector<std::string> toks;
                        Tokens::split(t[j], ",", toks);
                        
                        if (toks.size() == 1)
                        {
                            d.readR = stod(toks[0]);
                            d.readV = stod(toks[0]);
                        }
                        else
                        {
                            d.readR = stod(toks[0]);
                            d.readV = stod(toks[1]);
                        }
                    }
                    else if (formats[j] == "DP")
                    {
                        d.depth = stod(t[j]);
                    }
                }

                if (d.alt == ".")
                {
                    continue;
                }
                
                /*
                 * Parsing QUAL probability (TODO: Fix this...)
                 */
                
                //d.p = 1.0 - (exp(stold(fields[Field::QUAL]) / -10.0));
                
                f(d, p);
            }
        }
    };
}

#endif