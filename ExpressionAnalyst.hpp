#ifndef AS_EXPRESSION_ANALYST_HPP
#define AS_EXPRESSION_ANALYST_HPP

#include <string>
#include "types.hpp"

struct ExpressionStats
{
};

struct ExpressionAnalyst
{
	static ExpressionStats analyze(const std::string &file);
};

#endif