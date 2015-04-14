#include <map>
#include <iostream>
#include <assert.h>
#include "aligner.hpp"
#include "biology.hpp"
#include "expression.hpp"
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

    Feature f1, f2;
    
    /*
     * At the alignment, we're only interested in counting for the gene-level.
     * Isoform-level is another possibiltiy but there wouldn't be information
     * to distinguish ambiguous reads from alternative splicing.
     */
    
    std::map<GeneID, Counts> counts;

    std::for_each(r.seqs_gA.begin(), r.seqs_gA.end(), [&](const std::pair<GeneID, Sequins> &p)
    {
        counts[p.first] = 0;
    });
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if ((options.mode == ExonAlign   && align.spliced) ||
            (options.mode == SpliceAlign && !align.spliced))
        {
            return;
        }

		if (align.mapped)
		{
            const bool mapped = contains(r, align);

			if (align.id == r.id)
			{
				if (mapped)
				{
                    /*
                     * It's not enough that the read is mapped, it must also be mapped correctly.
                     * If the read maps to an exon, check whether it's mapped to an exon in the
                     * reference. Otherwise the read is spliced, check whether it maps to a
                     * junction in the reference.
                     */

                    const bool correct = (!align.spliced && find(r.fs.begin(), r.fs.end(), align, f1)) ||
                                          (align.spliced && checkSplice(r, align, f1, f2));

					if (correct)
					{
                        
                        assert(r.iso2Gene.count(f1.iID));
                        counts[r.iso2Gene.at(f1.iID)]++;
                        
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
				if (mapped)
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

    Expression::print(counts);    
    
    const auto cr = Expression::analyze(counts);
    
    // Either the samples are independent or the least detectable-abundant sequin is known
    assert(!cr.limit_count || r.seqs_gA.count(cr.limit_key));

    stats.sens.id = cr.limit_key;
    stats.sens.counts = cr.limit_count;
    stats.sens.abundance = cr.limit_count
                         ? r.seqs_gA.at(cr.limit_key).r.reads + r.seqs_gA.at(cr.limit_key).v.reads: NAN;

    if (options.writer)
    {
        options.writer->write(
            (boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%")
                % "diluation" % "tp" % "tn" % "fp" % "fn" % "sensitivity").str());
        options.writer->write(
            (boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%")
                % stats.dilution % stats.m.tp % stats.m.tn % stats.m.fp % stats.m.fn
                    % stats.sens.abundance).str());
    }

	return stats;
}