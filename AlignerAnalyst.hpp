#ifndef AS_ALIGNER_ANALYST_HPP
#define AS_ALIGNER_ANALYST_HPP

#include <string>
#include "types.hpp"

struct AlignerStats
{
	// Percentage of reads aligned with the silico
	Percentage p_si;

	// Percentage of reads aligned with the real samples
	Percentage p_sa;

	// Total number of reads aligned
	Reads n = 0;

	// Number of reads aligned to the silico chromosome
	Reads n_si = 0;

	// Number of reads aligned to the real sample
	Reads n_sa = 0;

	Percentage dilution;
};

struct AlignerAnalyst
{
	static AlignerStats analyze(const std::string &file);
};

#endif