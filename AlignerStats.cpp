#include <map>
#include <assert.h> 
#include "types.hpp"
#include "ParserFA.hpp"
#include "ParserSAM.hpp"
#include "AlignerCompare.hpp"

using namespace std;

AlignerStatistics AlignerCompare::analyze(const std::string &file)
{
	AlignerStatistics stats;

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
	 *
	 * In this context, positive refers to alignments mapped to the in-sillico chromosome, negative refers to alignments
	 * mapped but not to the in-sillico chromosome.
	 */

	ParserSAM::read(file, [&](const Alignment &align)
	{
		if (align.mapped)
		{
			// Mapped to the in-sillico chromosome
			if (align.id == sillico.id)
			{
				/*
				 * Is this a correct alignment? Assume it's for now...
				 */

				tp++;

				// if not -> fp++

				stats.n_si++;
			}

			// Not mapped to the chromosome (not necessarily mapped to the real-samples)
			else
			{
				// if really sequin then fn++
				// if really NOT sequin then tn++

				stats.n_sa++;
			}

			stats.n++;
		}
	});

	stats.p_si = static_cast<Percentage>(stats.n_si / stats.n);
	stats.p_sa = static_cast<Percentage>(stats.n_sa / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_si / stats.n_sa);

	return stats;
}