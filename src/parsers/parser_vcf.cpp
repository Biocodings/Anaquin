#include <map>
#include <assert.h>
#include "file.hpp"
#include "tokens.hpp"
#include "parsers/parser_vcf.hpp"
#include <iostream>
using namespace Spike;

/*
 * Refer to http://samtools.github.io/hts-specs/VCFv4.1.pdf for more details
 */

enum VCFField
{
    Chromo,
    Pos,
    ID,
    Ref,
    Alt,
    Qual,
    Filter,
    Info,
    Format,
    Format_Data,
};

static const std::map<std::string, Genotype> allele =
{
    { "0/0", HomozygousRef } , { "1/1", HomozygousAlt } , { "0/1", Heterzygous }
};

void ParserVCF::parse(const std::string &file, std::function<void (const VCFVariant &, const ParserProgress &)> fp)
{
    std::string line;
    File f(file);

    VCFVariant v;
    ParserProgress p;

    std::vector<std::string> t;
    std::vector<std::string> infos;
    std::vector<std::string> fields;
    std::vector<std::string> formats;

    while (f.nextLine(line))
    {
        p.i++;

		if (line[0] == '#')
		{
			continue;
		}

        Tokens::split(line, "\t", fields);

        v.id = fields[Chromo];
        v.l.start = v.l.end = stod(fields[Pos]);

        /*
         * Each base must be one of A,C,G,T,N (case insensitive). The value in the POS field refers
         * to the position of the first base in the string.
         */

        v.r = fields[Ref];

        /*
         * Additional information
         *
         *    AC: allele count in genotypes, for each ALT allele, in the same order as listed
         *    AF: allele frequency for each ALT allele in the same order as listed
         *    AN: total number of alleles in called genotypes
         *    DP: combined depth across samples
         */

        Tokens::split(fields[Info], ";", infos);

        for (const auto info : infos)
        {
            Tokens::split(info, "=", t);
            assert(t.size() == 2);

            if (t[0] == "AC") { v.ac = stof(t[1]); }
            if (t[0] == "AF") { v.af = stof(t[1]); }
            if (t[0] == "AN") { v.an = stof(t[1]); }
            if (t[0] == "DP") { v.dp = stod(t[1]); }
        }

        Tokens::split(fields[Format], ":", formats);

        /*
         * Comma separated list of alternate non-reference alleles called on at least one of the samples
         */

        std::vector<Sequence> alts;
        Tokens::split(fields[Alt], ",", alts);

        for (auto i = 0; i < alts.size(); i++)
        {
            v.a = alts[i];

            if (v.r.size() == v.a.size() && v.r.size() == 1)
            {
                v.m = SNP;
            }
            else
            {
                v.m = Indel;
            }

            Tokens::split(fields[Format_Data + i], ":", t);
            assert(t.size() == formats.size());

            for (auto j = 0; j < t.size(); j++)
            {
                if (formats[j] == "GT")
                {
                    v.gt = allele.at(t[j]);
                }
            }

            fp(v, p);
        }
	}
}