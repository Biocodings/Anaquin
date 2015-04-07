#include <vector>
#include <fstream>
#include <sstream>
#include "tokens.hpp"
#include "parser_fa.hpp"

using namespace Spike;

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
            
            Tokens::split(line, "|", tokens);

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
