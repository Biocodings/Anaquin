#ifndef AS_ALIGNER_ANALYST_HPP
#define AS_ALIGNER_ANALYST_HPP

#include <string>
#include "types.hpp"

struct AlignerStats
{
	// Percentage of reads aligned with the silico
	Percentage p_chromo;

	// Percentage of reads aligned with the real samples
	Percentage p_sample;

	// Total number of reads aligned
	Reads n = 0;

	// Number of reads aligned to the chromosome
	Reads n_chromo = 0;

	// Number of reads aligned to the real sample
	Reads n_sample = 0;

	Reads tp = 0; // Reads that are true-positive
	Reads tn = 0; // Reads that are true-negative
	Reads fp = 0; // Reads that are false-positive
	Reads fn = 0; // Reads that are false-negative

	Percentage sp; // Sensitivity
	Percentage sn; // Specificity

	Percentage dilution;
};

struct AlignerAnalyst
{
	static AlignerStats analyze(const std::string &file);
};

#endif