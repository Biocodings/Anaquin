#ifndef GI_PARSER_SAM_HPP
#define GI_PARSER_SAM_HPP

#include <functional>
#include "alignment.hpp"

struct ParserSAM
{
	static bool parse(const std::string &file, std::function<bool (const Alignment &)>);
};

#endif
