#ifndef AS_EXPRESSION_ANALYST_HPP
#define AS_EXPRESSION_ANALYST_HPP

#include "Types.hpp"
#include "Sequins.hpp"

struct ExpressionStats
{
};

struct ExpressionAnalyst
{
	static ExpressionStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif