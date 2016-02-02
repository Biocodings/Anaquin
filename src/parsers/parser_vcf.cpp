#include <map>
#include <iostream>
#include <assert.h>
#include "data/reader.hpp"
#include "data/tokens.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

/*
 * Refer to http://samtools.github.io/hts-specs/VCFv4.1.pdf for more details
 */

enum VCFField
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

static const std::map<std::string, Genotype> allele =
{
    { "0/0", HomozygousR },
    { "1/1", HomozygousA },
    { "0/1", Heterzygous },
    { "1/2", Heterzygous },
};

void ParserVCF::parse(const Reader &r, Callback c)
{
    std::string line;

    VCFVariant v;
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
        v.chrID = fields[VCFField::CHROM];

        // Eg: D_1_3_R
        v.id = fields[VCFField::ID];
        
        v.l.start = v.l.end = stod(fields[VCFField::POS]);

        // Reference allele
        v.ref = fields[VCFField::REF];

        /*
         * Additional information
         *
         *    AC: allele count in genotypes, for each ALT allele, in the same order as listed
         *    AF: allele frequency for each ALT allele in the same order as listed
         *    AN: total number of alleles in called genotypes
         *    DP: combined depth across samples
         */

        if (fields[VCFField::INFO] != ".")
        {
            Tokens::split(fields[VCFField::INFO], ";", infos);
            
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
        
        Tokens::split(fields[VCFField::FORMAT], ":", formats);
        assert(std::find(formats.begin(), formats.end(), "GT") != formats.end());

        /*
         * Comma separated list of alternate non-reference alleles called on at least one of the samples
         */

        std::vector<Sequence> alts;
        Tokens::split(fields[VCFField::ALT], ",", alts);

        for (auto i = 0; i < alts.size(); i++)
        {
            v.alt = alts[i];

            Tokens::split(fields[VCFField::FORMAT_DATA + i], ":", t);
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

            c(v, p);
        }
	}
}