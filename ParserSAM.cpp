#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "ParserSAM.hpp"
#include <boost/algorithm/string.hpp>

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

        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("\t"));

        // It's unmapped if the bit is 1
        align.mapped = !(stoi(tokens[1]) & (1 << 2));

        align.id  = tokens[2];
		align.seq = tokens[9];

        const auto len = static_cast<BasePair>(tokens[9].size());
        const auto start = stoi(tokens[3]);
        
        align.loc.update(start, start + len);
        assert(len == align.loc.length());
        
		if (!x(align))
        {
            break;
        }
	}
    
	return true;
}