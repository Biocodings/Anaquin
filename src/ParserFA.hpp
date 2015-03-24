#ifndef AS_READER_FA_HPP
#define AS_READER_FA_HPP

#include "Types.hpp"
#include <functional>
#include "Sequence.hpp"

struct ParserFA
{
	static bool parse(const std::string &file, std::function<void(const Sequence &)> x);
};

#endif
