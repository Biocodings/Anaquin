#include <map>
#include <assert.h> 
#include "types.hpp"
#include "ParserFA.hpp"
#include "ParserSAM.hpp"
#include "SillicoFactory.hpp"
#include "AlignerAnalyst.hpp"

using namespace std;

AlignerStats AlignerAnalyst::analyze(const std::string &file)
{
	AlignerStats stats;

	Reads tp = 0;
	Reads tn = 0;
	Reads fp = 0;
	Reads fn = 0;

	Sequence sillico;

	/*
	 * Extract the name of the in-sillico chromosome. Assume only a single chromosome in the file.
	 */

	ParserFA::parse("C:\\Sources\\QA\\data\\standards\\ChrT.5.10.fa", [&](const Sequence &s)
	{
		sillico = s;
	});

	/*
	 * Calculate the sensitivity and specificity for the experiment. Sensitivity and specificity are statistical measures
	 * of the performance of a binary test. Sensitivity (also known as true positive rate) measures the proportion of
	 * actual positives which are correctly identified as such (eg: the percentage of reads aligned to in-sillico correctly).
	 * Specificity (also known as true negative rate) measures the proportion of negatives which are correctly identified as
	 * such (eg: the percentage of reads fails to aligned to in-sillico while those reads come from the real-sample).
	 */

	ParserSAM::read(file, [&](const Alignment &align)
	{
		if (!align.mapped)
		{
			// A sequin fails to be mapped
			if (align.id == sillico.id)
			{
				// It's a false-negative because the mapping fails but it shouldn't
				fn++;
			}

			// 
			else
			{
				// It's a true-negative because the read doesn't belong to the sequins
				tn++;
			}

			stats.n_sa++;
		}
	});

	/*
	 * ????
	 */

	ParserFA::parse("C:\\Sources\\QA\\Data\\Standards\\RNAsequins.fa", [&](const Sequence &s)
	{
	});

	stats.p_si = static_cast<Percentage>(stats.n_si / stats.n);
	stats.p_sa = static_cast<Percentage>(stats.n_sa / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_si / stats.n_sa);

	return stats;
}