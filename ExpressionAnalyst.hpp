#ifndef AS_EXPRESSION_ANALYST_HPP
#define AS_EXPRESSION_ANALYST_HPP

#include "Types.hpp"
#include "Sequins.hpp"

struct ExpressionStats
{
    // Correlation for the samples
    double r;

    // Adjusted R2 for the linear model
    double r2;

    // Coefficient for the linear model
    double slope;
};

struct ExpressionAnalyst
{
	static ExpressionStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif