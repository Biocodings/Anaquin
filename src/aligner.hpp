#ifndef GI_ALIGNER_HPP
#define GI_ALIGNER_HPP

#include "sequins.hpp"
#include "confusion_matrix.hpp"

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

struct AlignerOptions
{
    AlignerOptions() {}
    AlignerOptions(Reads n) : n(n) {}

    // The number of reads to be analyzed
    Reads n = std::numeric_limits<Reads>::max();
};

struct Aligner
{
    // Analyze the aligner for the base-level
    static AlignerStats analyze(const std::string &file, AlignerOptions options = AlignerOptions());

    // Analyze for the aligner for the spliced-junction level
    static AlignerStats spliced(const std::string &file, AlignerOptions options = AlignerOptions());
};

#endif