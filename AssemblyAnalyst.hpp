#ifndef AS_ASSEMBLY_ANALYST_HPP
#define AS_ASSEMBLY_ANALYST_HPP

#include "ConfusionMatrix.hpp"

struct AssemblyStats
{
	/*
	 * Metrics for the base-level
	 */

	Percentage fp, tp, fn, tn;

	/*
	 * Metrics for exons
	 */

	Percentage e_fp, e_tp, e_fn, e_tn;
};

struct AssemblyAnalyst
{
	static AssemblyStats analyze(const std::string &file);
};

#endif