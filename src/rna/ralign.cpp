#include <iostream>
#include <assert.h>
#include "ralign.hpp"
#include "biology.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
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

RAlignStats RAlign::analyze(const std::string &file, const RAlign::Options &options)
{
    RAlignStats stats;
    const auto &r = Standard::instance();

    Feature f1, f2;

    /*
     * We're only interested in counting for the gene-level. Isoform-level is another possibiltiy
     * but there wouldn't be information to distinguish ambiguous reads from alternative splicing.
     */
    
    auto cb = countsForGenes();
    auto ce = countsForGenes();
    
    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if ((options.level == Exon   && align.spliced) ||
            (options.level == Splice && !align.spliced))
        {
            return;
        }

		if (align.mapped)
		{
            classify(stats, align, [&](const Alignment &)
            {
                if ((!align.spliced && find(r.fs.begin(), r.fs.end(), align, f1)) ||
                     (align.spliced && checkSplice(r, align, f1, f2)))
                {
                    assert(r.iso2Gene.count(f1.iID));
                    
                    ce[r.iso2Gene.at(f1.iID)]++;
                    cb[r.iso2Gene.at(f1.iID)]++;

                    return true;
                }
                else
                {
                    return false;
                }
            });
        }
    });

    assert(stats.nr + stats.nq == stats.n);
    
    const auto rb = Expression::analyze(cb);
    const auto re = Expression::analyze(ce);

    stats.sb = rb.sens(r.r_seqs_gA);
    stats.se = re.sens(r.r_seqs_gA);

    /*
     * Base-level statistics
     */

    AnalyzeReporter::reportSS("align.stats", stats, options.writer);

    /*
     * Exon-level statistics
     */

    /*
     * Counting statistics
     */
    
    AnalyzeReporter::reportCounts("base.counts", cb, options.writer);
    AnalyzeReporter::reportCounts("exon.counts", ce, options.writer);

	return stats;
}