#ifndef AS_ALIGNER_ANALYST_HPP
#define AS_ALIGNER_ANALYST_HPP

#include "Sequins.hpp"
#include "ConfusionMatrix.hpp"

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

	ConfusionMatrix m;

	Percentage dilution;
};

struct AlignerAnalyst
{
    static AlignerStats analyze(const std::string &file, Sequins s = Sequins(), Reads n = std::numeric_limits<Reads>::max());
};

#endif