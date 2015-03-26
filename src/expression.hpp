#ifndef GI_EXPRESSION_HPP
#define GI_EXPRESSION_HPP

#include "types.hpp"
#include "sequins.hpp"

struct ExpressionStats
{
    // Correlation for the samples
    double r;

    // Adjusted R2 for the linear model
    double r2;

    // Coefficient for the linear model
    double slope;
};

enum ExpressionMode
{
    GeneExpress,
    IsoformExpress,
};

struct Expression
{
	static ExpressionStats analyze(const std::string &file, ExpressionMode mode, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif