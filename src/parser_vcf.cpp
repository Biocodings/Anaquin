#include <map>
#include <vector>
#include <fstream>
#include "parser_vcf.hpp"
#include <boost/algorithm/string.hpp>

/*
 * Refer to http://samtools.github.io/hts-specs/VCFv4.1.pdf for more details
 */

enum VCFField
{
    CHROMO,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO,
    FORMAT,
    FORMAT_DATA_1
};

static std::map<std::string, AlleleType> alleleParser =
{
    { "0/0", HomozygousRef } , { "1/1", HomozygousAlt } , { "0/1", Heterzygous }
};

void ParserVCF::parse(const std::string &file, VCFHeaderF fh, VCFVariantF fv)
{
    std::string line;
    std::ifstream in(file);

    VCFVariant v;
    std::vector<std::string> t1;
    std::vector<std::string> t2;
    std::vector<std::string> t3;

    while (std::getline(in, line))
    {
		if (line.empty() || line[0] == '#')
		{
			continue;
		}

        boost::split(t1, line, boost::is_any_of("\t"));

        v.chromoID = t1[CHROMO];
        v.pos = stod(t1[POS]);
        v.varID = t1[ID];
        v.ref = t1[REF];

        /*
         * Example:
         *
         *     G,T
         */
        
        v.alts.clear();
        boost::split(v.alts, t1[ALT], boost::is_any_of(","));

        /*
         * Example:
         *
         *     GT:AD:DP:GQ:PL
         */
        
        t2.clear();
        boost::split(t2, t1[FORMAT], boost::is_any_of(":"));
        
        /*
         * Example:
         *
         *     1/1:0,128:128:99:5206,436,0
         */
        
        t3.clear();
        boost::split(t3, t1[FORMAT_DATA_1], boost::is_any_of(":"));

        assert(t2.size() == t3.size());
        
        for (auto i = 0; i < t2.size(); i++)
        {
            if (t2[i] == "GT")
            {
                v.type = alleleParser.at(t3[i]);
            }
        }
        
        fv(v);
	}
}