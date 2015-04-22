#include <map>
#include <assert.h>
#include "file.hpp"
#include "tokens.hpp"
#include "parsers/parser_vcf.hpp"

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
    Format_Data_1,
    Format_Data_2,
    Format_Data_3,
    Format_Data_4
};

void ParserVCF::parse(const std::string &file, VCFVariantF fv)
{
    std::string line;
    File f(file);

    VCFVariant v;

    std::vector<std::string> t1;
    std::vector<std::string> t2;
    std::vector<std::string> t3;

    while (f.nextLine(line))
    {
		if (line[0] == '#')
		{
			continue;
		}

        Tokens::split(line, "\t", t1);

        /*
         * An identifier from the reference genome. All entries for a specific CHROM should form a
         * contiguous block within the VCF file.
         */

        v.id = t1[Chromo];

        /*
         * The reference position, with the 1st base having position 1.
         */

        v.l.end = v.l.start = stod(t1[Pos]);

        /*
         * Each base must be one of A,C,G,T,N (case insensitive). The value in the POS field refers
         * to the position of the first base in the string.
         */
        
        v.r = t1[Ref];

        /*
         * Comma separated list of alternate non-reference alleles called on at least one of the samples
         */
        
        v.alts.clear();
        Tokens::split(t1[Alt], ",", v.alts);

        /*
         * Additional information
         */

        Tokens::split(t1[Info], ";", t2);
        
        for (const auto i : t2)
        {
            Tokens::split(i, "=", t3);

            if (t3.size() != 2)
            {
                throw std::runtime_error("Malformed VCF " + file);
            }

            v.info[t3[0]] = t3[1];
        }
        
        const std::map<std::string, Zygosity> allele =
        {
            { "0/0", Homozygous } , { "1/1", Homozygous } , { "0/1", Heterzygous }
        };

        /*
         * If genotype information is present, then the same types of data must be present for
         * all samples. First a FORMAT field is given specifying the data types and order.
         */

        Tokens::split(t1[Format], ":", t2);
        Tokens::split(t1[Format_Data_1], ":", t3);

        assert(t2.size() == t3.size());

        for (auto i = 0; i < t2.size(); i++)
        {
            if (t2[i] == "GT")
            {
                v.zy = allele.at(t3[i]);
            }
        }

        fv(v);
	}
}