#include <vector>
#include <fstream>
#include <sstream>
#include "types.hpp"
#include "ParserSAM.hpp"

using namespace std;

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

bool ParserSAM::read(const std::string &file, std::function<void(const Alignment &)> x)
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

        auto tokens = split(line, '\t');
        
        // It's unmapped if the bit is 1
        align.mapped = !(stoi(tokens[1]) & (1 << 2));
        
        align.id  = tokens[2];
		align.seq = tokens[9];
        align.pos = stoi(tokens[3]);

		x(align);
	}
    
	return true;
}