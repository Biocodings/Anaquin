#include <map>
#include <iostream>
#include <assert.h>
#include "aligner.hpp"
#include "biology.hpp"
#include "standard.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
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

AlignerStats Aligner::analyze(const std::string &file, const Aligner::Options &options)
{
    AlignerStats stats;
    const auto &r = Standard::instance();

    Feature f1, f2;

    /*
     * At the alignment, we're only interested in counting for the gene-level.
     * Isoform-level is another possibiltiy but there wouldn't be information
     * to distinguish ambiguous reads from alternative splicing.
     */
    
    INIT_COUNTER(c);
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if ((options.level == LevelExon   && align.spliced) ||
            (options.level == LevelSplice && !align.spliced))
        {
            return;
        }

		if (align.mapped)
		{
            const bool detected = contains(r, align);

			if (align.id == r.id)
			{
				if (detected)
				{
                    /*
                     * It's not enough that the read is mapped, it must also be mapped correctly.
                     * If the read maps to an exon, check whether it's mapped within a boundary.
                     * Otherwise the read is spliced, check whether it maps to a spliced-junction.
                     */

                    const bool correct = (!align.spliced && find(r.fs.begin(), r.fs.end(), align, f1)) ||
                                          (align.spliced && checkSplice(r, align, f1, f2));

					if (correct)
					{
                        assert(r.iso2Gene.count(f1.iID));
                        c[r.iso2Gene.at(f1.iID)]++;
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
				if (detected)
				{
					stats.m.fn++; // We can't find this!
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
    
    const auto cr = Expression::analyze(c);

    // Either the samples are independent or the least detectable-abundant sequin is known
    assert(!cr.limit_count || r.seqs_gA.count(cr.limit_key));

    stats.sens.id = cr.limit_key;
    stats.sens.counts = cr.limit_count;
    stats.sens.exp = cr.limit_count ? r.seqs_gA.at(cr.limit_key).r.raw +
                                      r.seqs_gA.at(cr.limit_key).v.raw: NAN;

    const std::string format = "%1%\t%2%\t%3%\t%4%";
    
    options.writer->open("align.stats");
    options.writer->write((boost::format(format) % "dl" % "sp" % "sn" % "ss").str());
    options.writer->write((boost::format(format) % stats.dilution()
                                                 % stats.m.sp()
                                                 % stats.m.sn()
                                                 % stats.sens.exp).str());
    options.writer->close();

	return stats;
}