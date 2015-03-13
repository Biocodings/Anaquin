#include <vector>
#include <fstream>
#include <sstream>
#include "sam_reader.hpp"

using namespace std;

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool SAMReader::read(const std::string &file, std::function<void(const Alignment &alignment)> x)
{
    std::string line;
    std::ifstream in(file);

    Alignment align;
    
    while (std::getline(in, line))
    {
        auto tokens = split(line, '\t');
        
        // It's unmapped if the bit is 1
        align.mapped = !(stoi(tokens[1]) & (1 << 2));
        
        align.name = tokens[2];

		x(align);
	}
    
	return true;
}