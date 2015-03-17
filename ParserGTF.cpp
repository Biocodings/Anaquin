#include <fstream>
#include <assert.h>
#include "ParserGTF.hpp"
#include <iostream>
#include <boost/algorithm/string.hpp>

using namespace std;

bool ParserGTF::parse(const std::string &file, std::function<void(const Feature &)> x)
{
    std::string line;
    std::ifstream in(file);

    Feature f;

	/*
	 * Fields must be tab-separated. Also, all but the final field in each feature line must contain a value;
	 * "empty" columns should be denoted with a '.'. Please refer to the online documentation for more details.
	 *
	 *    1. seqname
	 *    2. sourcre
	 *    3. feature
	 *    4. start
	 *    5. end
	 *    6. score
	 *    7. strand
	 *    8. frame
	 *    9. attribute
	 */
    
    std::vector<std::string> tokens;
    std::vector<std::string> options;
    std::vector<std::string> nameValue;

    while (std::getline(in, line))
    {
        boost::split(tokens, line, boost::is_any_of("\t"));

        f.chromo = tokens[0];
        f.l.update(stoi(tokens[3]), stoi(tokens[4]));

		if (tokens[2] == "exon")
		{
            f.type = Exon;
		}
		else if (tokens[2] == "CDS")
		{
			f.type = CDS;
		}
		else if (tokens[2] == "start_codon")
		{
			f.type = StartCodon;
		}
        
        /*
         * Eg: "gene_id "R_5_3_R"; transcript_id "R_5_3_R";"
         */
        
        boost::split(options, tokens[8], boost::is_any_of(";"));
        
        for (auto option : options)
        {
            if (!option.empty())
            {
                boost::trim(option);
                boost::split(nameValue, option, boost::is_any_of(" "));

                if (nameValue.size() == 2)
                {
                    f.options[nameValue[0]] = nameValue[1];
                }
            }
        }
        
        x(f);
	}
    
	return true;
}