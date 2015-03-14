#include <vector>
#include <fstream>
#include <sstream>
#include "ParserFA.hpp"

using namespace std;

bool ParserFA::parse(const std::string &file, std::function<void(const Sequence &)> x)
{
    std::string line;
    std::ifstream in(file);

    Sequence s;
    
    while (std::getline(in, line))
    {
		if (line[0] == '>')
		{
			if (!s.id.empty())
			{
				x(s);
			}

			auto tokens = split(line, '|');
			s.id = tokens[0];
			s.value.clear();
		}
		else
		{
			s.value += line;
		}
	}

	x(s);    
	return true;
}