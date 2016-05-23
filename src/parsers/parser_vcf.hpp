#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include <map>
#include "data/tokens.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "parsers/parser.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef CalledVariant Data;

        /*
         * Refer to http://samtools.github.io/hts-specs/VCFv4.1.pdf for more details
         */
        
        enum Field
        {
            CHROM,
            POS,
            ID,
            REF,
            ALT,
            QUAL,
            FILTER,
            INFO,
            FORMAT,
            FORMAT_DATA,
        };

        static void parse(const Reader &r, std::function<void (const Data &, const ParserProgress &)> f)
        {
            const std::map<std::string, Genotype> allele =
            {
                { "0/0", HomozygousR },
                { "1/1", HomozygousA },
                { "0/1", Heterzygous },
                { "1/2", Heterzygous },
            };

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
                
                if (line[0] == '#')
                {
                    continue;
                }
                
                Tokens::split(line, "\t", fields);
                
                // Eg: chrT
                d.cID = fields[Field::CHROM];
                
                // Eg: D_1_3_R
                d.id = fields[Field::ID];
                
                d.l.start = d.l.end = stod(fields[Field::POS]);
                
                // Reference allele
                d.ref = fields[Field::REF];
                
                /*
                 * Additional information
                 *
                 *    AF: allele frequency for each ALT allele in the same order as listed
                 */
                
                if (fields[Field::INFO] != ".")
                {
                    Tokens::split(fields[Field::INFO], ";", infos);
                    
                    for (const auto &info : infos)
                    {
                        Tokens::split(info, "=", t);
                        assert(t.size() == 2);
                        
                        if (t[0] == "AF") { d.af = stof(t[1]); }
                    }
                }
                
                Tokens::split(fields[Field::FORMAT], ":", formats);
                assert(std::find(formats.begin(), formats.end(), "GT") != formats.end());
                
                /*
                 * Comma separated list of alternate non-reference alleles called on at least one of the samples
                 */
                
                std::vector<Sequence> alts;
                Tokens::split(fields[Field::ALT], ",", alts);

                for (auto i = 0; i < alts.size(); i++)
                {
                    d.alt = alts[i];

                    Tokens::split(fields[Field::FORMAT_DATA + i], ":", t);
                    assert(t.size() == formats.size());
                    
                    for (auto j = 0; j < t.size(); j++)
                    {
                        if (formats[j] == "AD")
                        {
                            std::vector<std::string> toks;
                            
                            // Eg: 143,16 or 143
                            Tokens::split(t[j], ",", toks);
                            
                            assert(toks.size() == 1 || toks.size() == 2);
                            
                            if (toks.size() == 2)
                            {
                                d.readR = stod(toks[0]);
                                d.readV = stod(toks[1]);
                            }
                            else
                            {
                                d.readR = stod(toks[0]);
                                d.readV = stod(toks[0]);
                            }
                        }
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