#ifndef AS_ASSEMBLY_COMPARE_HPP
#define AS_ASSEMBLY_COMPARE_HPP

#include <string>
#include "types.hpp"

struct AssemblyStatistics
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

struct AssemblyCompare
{
	static AssemblyStatistics analyze(const std::string &file);
};

#endif