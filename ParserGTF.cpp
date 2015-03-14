#include <vector>
#include <fstream>
#include <sstream>
#include "types.hpp"
#include "ParserGTF.hpp"

using namespace std;

bool ParserGTF::parse(const std::string &file, FeatureReader &reader)
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

    while (std::getline(in, line))
    {
        const auto tokens = split(line, '\t');        
        f.pos = stoi(tokens[3]);
        
		if (tokens[2] == "exon")
		{
			reader.exon(f);
		}

        reader.all(f);
	}
    
	return true;
}