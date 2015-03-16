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
        if (++i > n)
        {
            return false;
        }
        else if (align.id == r.id)
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
						stats.m.tp++;
					}
					else
					{
						stats.m.fp++;
					}
				}
				else
				{
					stats.m.fp++;
				}
			}
			else
			{
				if (matchChromoBoundary(r, align))
				{
					stats.m.fn++;
				}
				else
				{
					stats.m.tn++;
				}
			}
		}

        stats.n++;        
        return true;
    });
    
    assert(stats.n);
    
    std::cout << "Total reads: " << stats.n << std::endl;
    std::cout << "Reference reads: " << stats.nr << std::endl;
    std::cout << "Query reads: " << stats.nq << std::endl;
    
	std::cout << "TP: " << stats.m.tp << std::endl;
	std::cout << "FN: " << stats.m.fn << std::endl;
	std::cout << "FP: " << stats.m.fp << std::endl;

	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
    stats.dilution = stats.nq ? static_cast<Percentage>(stats.nr / stats.nq) : 1;

	return stats;
}