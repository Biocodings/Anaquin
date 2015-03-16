#ifndef AS_PARSER_SAM_HPP
#define AS_PARSER_SAM_HPP

#include <functional>
#include "alignment.hpp"

struct ParserSAM
{
	static bool read(const std::string &file, std::function<bool (const Alignment &)> x);
};

#endif