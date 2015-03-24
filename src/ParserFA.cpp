#include <vector>
#include <fstream>
#include <sstream>
#include "ParserFA.hpp"
#include <boost/algorithm/string.hpp>

bool ParserFA::parse(const std::string &file, std::function<void(const Sequence &)> x)
{
    std::string line;
    std::ifstream in(file);

    Sequence s;
    std::vector<std::string> tokens;
    
    while (std::getline(in, line))
    {
		if (line[0] == '>')
		{
			if (!s.id.empty())
			{
				x(s);
			}
            
            boost::split(tokens, line, boost::is_any_of("|"));

			// Remove the '<' prefix			
			s.id = tokens[0].substr(1, tokens[0].size());

			// Reset for the next sequence
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
