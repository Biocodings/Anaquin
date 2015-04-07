#include <iostream>
#include <assert.h>
#include "aligner.hpp"
#include "biology.hpp"
#include "parser_bed.hpp"
#include "parser_sam.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "standard_factory.hpp"

using namespace Spike;

AlignerStats Aligner::spliced(const std::string &file, const AlignerOptions &options)
{
    const auto r = StandardFactory::reference();

    AlignerStats stats;
    
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
        if (align.id == r.id)
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
                if (contains(r, align))
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
                //if (matchChromoBoundary(r, align))
                {
                    stats.m.fn++;
                }
                //else
                {
                    stats.m.tn++;
                }
            }

            stats.n++;
        }
    });

    return stats;
}

AlignerStats Aligner::analyze(const std::string &file, const AlignerOptions &options)
{
    AlignerStats stats;
    const auto r = StandardFactory::reference();

    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if (!options.spliced && align.spliced)
        {
            return;
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
				if (contains(r, align))
				{
					if (contains(r.fs.begin(), r.fs.end(), align))
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
				if (contains(r, align))
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
    });

    assert(stats.nr + stats.nq == stats.n);
    
	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
    
    // Proportion of the reads from the reference
    stats.dilution = stats.nq ? static_cast<Percentage>(stats.nr / stats.nq) : 1;

    Writer w("align.csv");
 
    w.write((boost::format("%1%\t%2%\t%3%\t%4%\t%5%") % "diluation" % "tp" % "tn" % "fp" % "fn").str());
    w.write((boost::format("%1%\t%2%\t%3%\t%4%\t%5%") % stats.dilution % stats.m.tp % stats.m.tn % stats.m.fp % stats.m.fn).str());
    
	return stats;
}