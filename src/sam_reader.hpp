#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include <functional>
#include "alignment.hpp"

struct SAMReader
{
	static bool read(const std::string &file, std::function<void(const Alignment &alignment)> x);
};

#endif