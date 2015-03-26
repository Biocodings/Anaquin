#include <vector>
#include <fstream>
#include <sstream>
#include "parser_fa.hpp"
#include <boost/algorithm/string.hpp>

void ParserFA::parse(const std::string &file, std::function<void (const FASequence &)> f)
{
    std::string line;
    std::ifstream in(file);

    FASequence s;
    std::vector<std::string> tokens;
    
    while (std::getline(in, line))
    {
		if (line[0] == '>')
		{
			if (!s.id.empty())
			{
				f(s);
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

	f(s);
}
