#include <vector>
#include <fstream>
#include <sstream>
#include "reader_fa.hpp"

using namespace std;

static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    
    return elems;
}

static std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool ReaderFA::read(const std::string &file, std::function<void(const Sequence &s)> x)
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
    
	return true;
}