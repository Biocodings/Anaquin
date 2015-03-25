#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "parser_vcf.hpp"
#include <boost/algorithm/string.hpp>

void ParserVCF::parse(const std::string &file, VCFHeaderF &fh, VCFVariantF &fv)
{
    std::string line;
    std::ifstream in(file);

    VCFVariant v;
    std::vector<std::string> tokens;

    while (std::getline(in, line))
    {
		if (line.empty() || line[0] == '#')
		{
			continue;
		}

        boost::split(tokens, line, boost::is_any_of("\t"));

        v.chromoID = tokens[0];
        v.pos = stod(tokens[1]);
        v.varID = tokens[2];
        v.ref = tokens[3];

        
        //Comma separated list of alternate non-reference alleles
        
        //boost::split(tokens, tokens[4], boost::is_any_of(","));

        fv(v);
	}
}