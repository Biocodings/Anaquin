#include <map>
#include <math.h>
#include <iostream>
#include <assert.h> 
#include "ParserSAM.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

static bool matchGeneBoundary(const Chromosome &r, const Alignment &align)
{
    for (auto f: r.fs)
    {
        if (align.start >= r.start && align.end <= r.end)
        {
            return true;
        }
    }
    
    return false;
}

static bool matchChromoBoundary(const Chromosome &r, const Alignment &align)
{
	return (align.start >= r.start && align.end <= r.end);
}

AlignerStats AlignerAnalyst::analyze(const std::string &file, Sequins s, Reads n)
{
    AlignerStats stats;

    // The reference silico chromosome
    const auto r = StandardFactory::reference();

    /*
     * Calculate the sensitivity and specificity for the experiment. They are statistical measures of the performance
     * of a binary test. Sensitivity (also known true-positive rate) measures the proportion of actual positives which
     * are correctly identified as such (eg: the percentage of reads aligned to the reference correctly). Specificity
     * (also known as true-negative rate) measures the proportion of negatives which are correctly identified as such
     * (eg: the percentage of reads fails to align to the reference incorrectly).
     */
    
    Reads i = 0;
    
    ParserSAM::read(file, [&](const Alignment &align)
    {
        if (i++ > n)
        {
            return false;
        }
        
        if (align.id == r.id)
        {
            stats.nr++;
        }
        else
        {
            stats.nq++;
        }
        
		if (align.mapped)
		{
			if (align.id == r.id)
			{
				if (matchChromoBoundary(r, align))
				{
					if (matchGeneBoundary(r, align))
					{
						stats.tp++;
					}
					else
					{
						stats.fp++;
					}
				}
				else
				{
					stats.fp++;
				}
			}
			else
			{
				if (matchChromoBoundary(r, align))
				{
					stats.fn++;
				}
				else
				{
					stats.tn++;
				}
			}
		}

        stats.n++;        
        return true;
    });
    
	stats.sp = (stats.tp + stats.fn) ? stats.tp / (stats.tp + stats.fn) : NAN;
	stats.sn = (stats.fp + stats.tn) ? stats.tn / (stats.fp + stats.tn) : NAN;

    std::cout << "Total reads: " << stats.n << std::endl;
    std::cout << "Reference reads: " << stats.nr << std::endl;
    std::cout << "Query reads: " << stats.nq << std::endl;
    
    std::cout << "TP: " << stats.tp << std::endl;
    std::cout << "FN: " << stats.fn << std::endl;
    std::cout << "FP: " << stats.fp << std::endl;

	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
	stats.dilution = static_cast<Percentage>(stats.nr / stats.nq);

	return stats;
}