#ifndef AS_PARSER_SAM_HPP
#define AS_PARSER_SAM_HPP

#include <functional>
#include "Alignment.hpp"

struct ParserSAM
{
	static bool parse(const std::string &file, std::function<bool (const Alignment &)>);
};

#endif