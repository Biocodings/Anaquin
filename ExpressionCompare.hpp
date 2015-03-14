#ifndef AS_EXPRESSION_COMPARE_HPP
#define AS_EXPRESSION_COMPARE_HPP

#include <string>
#include "types.hpp"

struct ExpressionStatistics
{
};

struct ExpressionCompare
{
	static ExpressionStatistics analyze(const std::string &file);
};

#endif