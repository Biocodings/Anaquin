#include <map>
#include <assert.h> 
#include "types.hpp"
#include "ParserFA.hpp"
#include "ParserSAM.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

AlignerStats AlignerAnalyst::analyze(const std::string &file)
{
	AlignerStats stats;

	// Nmae of the chromosome, we'll nee the name for comparison
	const auto chromo = StandardFactory::chromoName();

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
			if (align.id == chromo)
			{
				stats.n_chromo++;
                
				// It's a false-negative because the mapping fails but it shouldn't
				stats.fn++;
			}
			else
			{
				stats.n_sample++;
                
				// It's a true-negative because the read doesn't belong to the sequins
				stats.tn++;
			}
		}
        else
        {
			if (align.id == chromo)
            {
				stats.n_chromo++;
                
                if (features.count(align.pos))
                {
                    // True-positive because the alignment to chromosome is correct
					stats.tp++;
                }
                else
                {
					stats.n_sample++;
                    
                    // Negative-positive because the alignment to chromosome is incorrect
					stats.fp++;
                }
            }
			else
			{
				stats.n_sample++;
			}
        }
        
        stats.n++;
	});

	stats.sp = (stats.tp + stats.fn) ? static_cast<Percentage>(stats.tp / (stats.tp + stats.fn)) : NAN;
	stats.sn = (stats.fp + stats.tn) ? static_cast<Percentage>(stats.tn / (stats.fp + stats.tn)) : NAN;

	stats.p_chromo = static_cast<Percentage>(stats.n_chromo / stats.n);
	stats.p_sample = static_cast<Percentage>(stats.n_sample / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_chromo / stats.n_sample);

	return stats;
}