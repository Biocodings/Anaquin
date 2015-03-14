#ifndef AS_READER_FA_HPP
#define AS_READER_FA_HPP

#include "types.hpp"
#include <functional>
#include "sequence.hpp"

struct ParserFA
{
	static bool parse(const std::string &file, std::function<void(const Sequence &)> x);
};

#endif