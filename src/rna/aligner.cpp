#include <map>
#include <iostream>
#include <assert.h>
#include "aligner.hpp"
#include "biology.hpp"
#include <ss/stats.hpp>
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "standard_factory.hpp"
#include "parsers/parser_bed.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

static bool checkSplice(const Standard &r, const Alignment &align, Feature &e1, Feature &e2)
{
    assert(align.spliced);

    /*
     * A spliced alignment is correct if it aligns to two consecutive exons
     */
        
    for (auto i = 0; i < r.exons.size(); i++)
    {
        const auto &exon_1 = r.exons[i];
        
        // Check if it starts inside exon_1
        if (align.l.start >= exon_1.l.start && align.l.start <= exon_1.l.end && align.l.end > exon_1.l.end
            && i != r.exons.size() - 1)
        {
            for (auto j = i + 1; j < r.exons.size(); j++)
            {
                const auto exon_2 = r.exons[j];
                
                if (exon_1.iID == exon_2.iID)
                {
                    // Check if it ends inside exon_2
                    if (align.l.start < exon_2.l.start && align.l.end >= exon_2.l.start && align.l.end <= exon_2.l.end)
                    {
                        e1 = exon_1;
                        e2 = exon_2;
                        return true;
                    }
                }
            }
        }
    }
    
    return false;
}

AlignerStats Aligner::analyze(const std::string &file, const AlignerOptions &options)
{
    AlignerStats stats;
    const auto r = StandardFactory::reference();

    // Used for alternative splicing
    Feature e1, e2;
    
    std::map<TranscriptID, unsigned> counts;
    std::for_each(r.seqs_iA.begin(), r.seqs_iA.end(), [&](const std::pair<TranscriptID, Sequin> &p)
    {
        counts[p.first] = 0;
    });
    
    // Make sure we have an entry for each protein isoform
    assert(counts.size() == r.seqs_iA.size());

    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if ((options.mode == ExonAlign   && align.spliced) ||
            (options.mode == SpliceAlign && !align.spliced))
        {
            return;
        }

		if (align.mapped)
		{
            // Whether the read has mapped to the reference (by checking via the locus)
            const bool ref_mapped = contains(r, align);

			if (align.id == r.id)
			{
				if (ref_mapped)
				{
                    const bool detected = !align.spliced || (align.spliced && checkSplice(r, align, e1, e2));

					if (detected)
					{
                        Feature f;
                        
                        if (!align.spliced && !find(r.fs.begin(), r.fs.end(), align, f))
                        {
                            throw std::runtime_error("Not splicing, not found????");
                        }
                                                
                        if (align.spliced)
                        {
                            f = e1;
                        }
                        
                        if (!counts.count(f.iID))
                        {
                            if (f.iID == "R_5_3_V")
                            {
                                assert(counts.count("R_5_3_R"));
                                counts["R_5_3_R"]++;
                            }
                            else
                            {
                                throw std::runtime_error("Not splicing, not found????");
                            }
                        }
                        else
                        {
                            assert(counts.count(f.iID));
                            counts[f.iID]++;
                        }
                        
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
				if (ref_mapped)
				{
					stats.m.fn++;
				}
				else
				{
					stats.m.tn++;
				}
			}
		}

        if (align.id == r.id)
        {
            stats.nr++;
        }
        else
        {
            stats.nq++;
        }
        
        stats.n++;        
    });

    assert(stats.nr + stats.nq == stats.n);
    
	stats.pr = static_cast<Percentage>(stats.nr / stats.n);
	stats.pq = static_cast<Percentage>(stats.nq / stats.n);
    
    // Proportion of the reads from the reference
    stats.dilution = stats.nq ? static_cast<Percentage>(stats.nr / stats.nq) : 1;

    /*
     * The counts for each sequin is needed to calculate the limit of sensitivity.
     */

    const auto cr = SS::analyze(counts);
    assert(r.seqs_iA.count(cr.min_k));
    
    // The least abdunance while still detectable in this experiment
    const auto limit_sensitivity = r.seqs_iA.at(cr.min_k).reads;
    
    if (options.writer)
    {
        options.writer->write(
            (boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%")
                % "diluation" % "tp" % "tn" % "fp" % "fn" % "detection").str());
        options.writer->write(
            (boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%")
                % stats.dilution % stats.m.tp % stats.m.tn % stats.m.fp % stats.m.fn
                    % limit_sensitivity).str());
    }

	return stats;
}