#include <math.h>
#include <iostream>
#include <assert.h>
#include <limits>
#include "aligner.hpp"
#include "biology.hpp"
#include "parser_bed.hpp"
#include "parser_sam.hpp"
#include "statistics.hpp"
#include "standard_factory.hpp"

static bool matchGeneBoundary(const Standard &r, const Alignment &align)
{
    for (auto f: r.fs) // TODO: Double check no intron is here.... // TODO: Check spliced reads
    {
        if (r.l.contains(align.l))
        {
            return true;
        }
    }
    
    return false;
}

static bool matchChromoBoundary(const Standard &r, const Alignment &align)
{
    return r.l.contains(align.l);
}

AlignerStats Aligner::spliced(const std::string &file, AlignerOptions options)
{
    const auto r = StandardFactory::reference();

    AlignerStats stats;
    
    Reads i = 0;
    
    // Check whether a spliced-alignment is correct
    auto f = [&](const Alignment &align)
    {
        /*
         * An spliced alignment is correct if it aligns to two consecutive exons
         */
     
        for (auto i = 0; i < r.exons.size(); i++)
        {
            const auto &exon_1 = r.exons[i];

            // Check if it starts inside exon_1
            if (align.l.start >= exon_1.l.start && align.l.start < exon_1.l.end && align.l.end > exon_1.l.end
                    && i != r.exons.size() - 1)
            {
                // Any ... ?
                
                for (auto j = i + 1; j < r.exons.size(); j++)
                {
                    const auto exon_2 = r.exons[j];
                    
                    if (exon_1.iID == exon_2.iID)
                    {
                        // Check if it ends inside exon_2
                        if (align.l.start < exon_2.l.start && align.l.end > exon_2.l.start && align.l.end < exon_2.l.end)
                        {
                            return true;
                        }
                    }
                }
            }
        }
        
        return false;
    };

    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if (++i > options.n)
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
        
        if (align.mapped && align.spliced)
        {
            if (align.id == r.id)
            {
                if (contains__(r, align))
                {
                    if (f(align))
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

            stats.n++;
        }
        
        return true;
    });

    return stats;
}

AlignerStats Aligner::analyze(const std::string &file, AlignerOptions options)
{
    AlignerStats stats;
    const auto r = StandardFactory::reference();

    /*
     * Calculate the sensitivity and specificity for the experiment. They are statistical measures of the performance
     * of a binary test. Sensitivity (also known true-positive rate) measures the proportion of actual positives which
     * are correctly identified as such (eg: the percentage of reads aligned to the reference correctly). Specificity
     * (also known as true-negative rate) measures the proportion of negatives which are correctly identified as such
     * (eg: the percentage of reads fails to align to the reference incorrectly).
     */
    
    Reads i = 0;
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if (++i > options.n)
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
                    // TODO: Rename to feature
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
    
	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
    stats.dilution = stats.nq ? static_cast<Percentage>(stats.nr / stats.nq) : 1;

    std::cout << "Dilution: " << stats.dilution << std::endl;
    std::cout << "TP: " << stats.m.tp << std::endl;
    std::cout << "TN: " << stats.m.tn << std::endl;
    std::cout << "FP: " << stats.m.fp << std::endl;
    std::cout << "FN: " << stats.m.fn << std::endl;
    
	return stats;
}
