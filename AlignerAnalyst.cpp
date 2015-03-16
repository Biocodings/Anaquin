#include <map>
#include <math.h>
#include <assert.h> 
#include "ParserSAM.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

static bool matchAlignWithRef(const Chromosome &r, const Alignment &align)
{
    for (auto f: r.fs)
    {
        if (f.start == align.start)
        {
            return true;
        }
    }
    
    return false;
}

AlignerStats AlignerAnalyst::analyze(const std::string &file)
{
	AlignerStats stats;

    // The reference chromosome
	const auto r = StandardFactory::reference();

    /*
     * Calculate the sensitivity and specificity for the experiment. They are statistical measures of the performance
     * of a binary test. Sensitivity (also known true-positive rate) measures the proportion of actual positives which
     * are correctly identified as such (eg: the percentage of reads aligned to the reference correctly). Specificity
     * (also known as true-negative rate) measures the proportion of negatives which are correctly identified as such
     * (eg: the percentage of reads fails to align to the reference incorrectly).
     *
     * The calculations is very similar to Cuffcompare in the Cufflinks package. Refer to reportStats() in cuffcompare.cpp
     * for more details.
     */
    
    ParserSAM::read(file, [&](const Alignment &align)
    {
        if (align.id == r.id)
        {
            stats.n_r++;
        }
        else
        {
            stats.n_q++;
        }
        
        if (align.id == r.id)
        {
            if (align.mapped)
            {
                if (matchAlignWithRef(r, align))
                {
                    stats.tp++;
                }
                else
                {
                    stats.fn++;
                }
            }
            else
            {
                stats.fn++;
            }
        }
        
        stats.n++;
    });
    
    stats.fp = stats.n_r - stats.tp;
    stats.sp = (stats.tp + stats.fp) ? stats.tp / (stats.tp + stats.fp) : NAN;
    stats.sn = (stats.tp + stats.fn) ? stats.tp / (stats.tp + stats.fn) : NAN;

	stats.p_r = static_cast<Percentage>(stats.n_r / stats.n);
	stats.p_q = static_cast<Percentage>(stats.n_q / stats.n);
	stats.dilution = static_cast<Percentage>(stats.n_r / stats.n_q);

	return stats;
}