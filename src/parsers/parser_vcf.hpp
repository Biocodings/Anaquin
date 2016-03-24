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

        static void parse(const Reader &r, std::function<void (const Data &, const ParserProgress &)> c)
        {
            const std::map<std::string, Genotype> allele =
            {
                { "0/0", HomozygousR },
                { "1/1", HomozygousA },
                { "0/1", Heterzygous },
                { "1/2", Heterzygous },
            };

            std::string line;
            
            Data v;
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
                v.chrID = fields[Field::CHROM];
                
                // Eg: D_1_3_R
                v.id = fields[Field::ID];
                
                v.l.start = v.l.end = stod(fields[Field::POS]);
                
                // Reference allele
                v.ref = fields[Field::REF];
                
                /*
                 * Additional information
                 *
                 *    AC: allele count in genotypes, for each ALT allele, in the same order as listed
                 *    AF: allele frequency for each ALT allele in the same order as listed
                 *    AN: total number of alleles in called genotypes
                 *    DP: combined depth across samples
                 */
                
                if (fields[Field::INFO] != ".")
                {
                    Tokens::split(fields[Field::INFO], ";", infos);
                    
                    for (const auto &info : infos)
                    {
                        Tokens::split(info, "=", t);
                        assert(t.size() == 2);
                        
                        if (t[0] == "AC") { v.ac = stof(t[1]); }
                        if (t[0] == "AF") { v.af = stof(t[1]); }
                        if (t[0] == "AN") { v.an = stof(t[1]); }
                        if (t[0] == "DP") { v.dp = stod(t[1]); }
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
                    v.alt = alts[i];

                    Tokens::split(fields[Field::FORMAT_DATA + i], ":", t);
                    assert(t.size() == formats.size());
                    
                    for (auto j = 0; j < t.size(); j++)
                    {
                        if (formats[j] == "GT")
                        {
                            v.gt = allele.at(t[j]);
                        }
                        else if (formats[j] == "AD")
                        {
                            std::vector<std::string> tokens;
                            
                            // Eg: 143,16 or 143
                            Tokens::split(t[j], ",", tokens);
                            
                            assert(tokens.size() == 1 || tokens.size() == 2);
                            
                            if (tokens.size() == 2)
                            {
                                v.dp_r = stod(tokens[0]);
                                v.dp_a = stod(tokens[1]);
                            }
                            else
                            {
                                v.dp_r = stod(tokens[0]);
                                v.dp_a = stod(tokens[0]);
                            }
                        }
                    }
                }

                c(v, p);
            }
        }
    };
}

#endif