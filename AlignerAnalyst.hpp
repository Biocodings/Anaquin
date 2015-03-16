#ifndef AS_ALIGNER_ANALYST_HPP
#define AS_ALIGNER_ANALYST_HPP

#include "types.hpp"

struct AlignerStats
{
	// Percentage of reads aligned with the reference chromosome
	Percentage p_r;

	// Percentage of reads aligned with the query samples
	Percentage p_q;

	// Total number of reads aligned
	Reads n = 0;

	// Number of reads aligned to the chromosome
	Reads n_r = 0;

	// Number of reads aligned to the real sample
	Reads n_q = 0;

    // Number of perfectly aligned reads (true positives)
	Reads tp = 0;

    // Number of reads that have been aligned incorrectly
	Reads fp = 0;

    // Number of reads that are in the chromosome but failed to align
    Reads fn = 0;

    // Sensitivity, the ability of the experiment to detect positively
	Percentage sp;
    
    // Specificity, the ability of the experiement to detect negatively
	Percentage sn;

	Percentage dilution;
};

struct AlignerAnalyst
{
	static AlignerStats analyze(const std::string &file);
};

#endif