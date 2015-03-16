#ifndef AS_ALIGNER_ANALYST_HPP
#define AS_ALIGNER_ANALYST_HPP

#include "Types.hpp"
#include "Sequins.hpp"

struct AlignerStats
{
	// Percentage of reads aligned with the reference chromosome
	Percentage pr;

	// Percentage of reads aligned with the query samples
	Percentage pq;

	// Total number of reads aligned
	Reads n = 0;

	// Number of reads aligned to the chromosome
	Reads nr = 0;

	// Number of reads aligned to the real sample
	Reads nq = 0;

    // Number of perfectly aligned reads (true-positives)
	Reads tp = 0;

	// Number of incorrectly aligned reads (false-positive)
	Reads fp = 0;

    // Number of reads that are not aligned to the chromosome correctly (true-negative)
	Reads tn = 0;

    // Number of reads that are not aligned to the chromosome incorrectly (false-negative)
    Reads fn = 0;

    // Sensitivity, the ability of the experiment to detect positively
	Percentage sp;
    
    // Specificity, the ability of the experiement to detect negatively
	Percentage sn;

	Percentage dilution;
};

struct AlignerAnalyst
{
    static AlignerStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif