#ifndef READER_FA_HPP
#define READER_FA_HPP

#include "types.hpp"
#include <functional>
#include "sequence.hpp"

struct ReaderFA
{
	static bool read(const std::string &file, std::function<void(const Sequence &s)> x);
};

#endif