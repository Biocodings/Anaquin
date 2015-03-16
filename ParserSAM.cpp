#include <vector>
#include <fstream>
#include <sstream>
#include "types.hpp"
#include "ParserSAM.hpp"

using namespace std;

bool ParserSAM::read(const std::string &file, std::function<bool (const Alignment &)> x)
{
    std::string line;
    std::ifstream in(file);

    Alignment align;
    
    while (std::getline(in, line))
    {
		if (line.empty() || line[0] == '@')
		{
			continue;
		}

        const auto tokens = split(line, '\t');

        // It's unmapped if the bit is 1
        align.mapped = !(stoi(tokens[1]) & (1 << 2));

        align.id  = tokens[2];
		align.seq = tokens[9];
        align.start = stoi(tokens[3]);
		align.length = static_cast<BasePair>(tokens[9].size());
		align.end = align.start + align.length;

		if (!x(align))
        {
            break;
        }
	}
    
	return true;
}