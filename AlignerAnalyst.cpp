#include <map>
#include <math.h>
#include <iostream>
#include <assert.h> 
#include "ParserSAM.hpp"
#include "AlignerAnalyst.hpp"
#include "StandardFactory.hpp"

using namespace std;

/*
 * Reference from Cuffcompare.cpp:
 *
 *   bool exon_match(GXSeg& r, GXSeg& q, uint fuzz=0) {
 *       uint sd = (r.start>q.start) ? r.start-q.start : q.start-r.start;
 *       uint ed = (r.end>q.end) ? r.end-q.end : q.end-r.end;
 *       uint ex_range=exonEndRange;
 *       if (ex_range<=fuzz) ex_range=fuzz;
 *       if ((r.flags&1) && (q.flags&1)) {
 * 	         if (sd>ex_range) return false;
 *       }
 *       else {
 * 	         if (sd>fuzz) return false;
 *       }
 *       if ((r.flags&2) && (q.flags&2)) {
 * 	         if (ed>ex_range) return false;
 *       }
 *       else {
 * 	         if (ed>fuzz) return false;
 *       }
 *       return true;
 *   }
 */

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

AlignerStats AlignerAnalyst::analyze(const std::string &file, Sequins s, Reads n)
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
                    stats.tp = stats.tp;
                }
            }
            else
            {
                stats.fn++;
            }
        }
        
        stats.n++;        
        return true;
    });
    
    stats.fp = stats.nr - stats.tp;
    stats.sp = (stats.tp + stats.fp) ? stats.tp / (stats.tp + stats.fp) : NAN;
    stats.sn = (stats.tp + stats.fn) ? stats.tp / (stats.tp + stats.fn) : NAN;

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