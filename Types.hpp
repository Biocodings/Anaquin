#ifndef AS_TYPES_HPP
#define AS_TYPES_HPP

#include <string>
#include <vector>
#include <sstream>

typedef std::string ID;
typedef long long Reads;
typedef float Percentage;

// Eg: 388488 from the first matching base
typedef long long Position;

inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}

	return elems;
}

inline std::vector<std::string> split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

#endif