#include <map>
#include <assert.h> 
#include "types.hpp"
#include "ParserFA.hpp"
#include "ParserSAM.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

#include <iostream>

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

	ParserFA::parse("/Users/user1/Sources/ABCD/standards/ChrT.5.10.fa", [&](const Sequence &s)
	{
		sillico = s;
	});
    
	const auto features = StandardFactory::features();

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
                stats.n_si++;
                
				// It's a false-negative because the mapping fails but it shouldn't
				fn++;
			}
			else
			{
                stats.n_sa++;
                
				// It's a true-negative because the read doesn't belong to the sequins
				tn++;
			}
		}
        else
        {
            if (align.id == sillico.id)
            {
                stats.n_si++;
                
                if (features.count(align.pos))
                {
                    // True-positive because the alignment to chromosome is correct
                    tp++;
                }
                else
                {
                    stats.n_sa++;
                    
                    // Negative-positive because the alignment to chromosome is incorrect
                    fp++;
                }
            }
        }
        
        stats.n++;
	});
    
	stats.p_si = static_cast<Percentage>(stats.n_si / stats.n);
	stats.p_sa = static_cast<Percentage>(stats.n_sa / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_si / stats.n_sa);

	return stats;
}