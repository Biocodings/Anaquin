#include <vector>
#include <fstream>
#include <sstream>
#include "Types.hpp"
#include "ParserGTF.hpp"

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

    while (std::getline(in, line))
    {
        const auto tokens = split(line, '\t');        

        f.id    = tokens[0];
        f.start = stoi(tokens[3]);
        f.end   = stoi(tokens[4]);

		if (tokens[2] == "exon")
		{
            f.type = Exon;
		}

        x(f);
	}
    
	return true;
}