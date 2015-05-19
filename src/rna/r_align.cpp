#include <assert.h>
#include "r_align.hpp"
#include "expression.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

static bool checkSplice(const Alignment &align, Feature &f)
{
    assert(align.spliced);
    const auto &s = Standard::instance();
    
    for (auto i = 0; i < s.r_introns.size(); i++)
    {
        if (align.l == s.r_introns[i].l)
        {
            f = s.r_introns[i];
            return true;
        }
    }

    return false;
}

RAlignStats RAlign::analyze(const std::string &file, const Options &options)
{
    RAlignStats stats;
    const auto &s = Standard::instance();

    std::vector<Alignment> q_exons, q_introns;

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        Feature f;
        
        if (align.id != s.id || !align.mapped)
        {
            return;
        }

        /*
         * Classify at the exon level
         */

        if (!align.spliced)
        {
            q_exons.push_back(align);

            if (classify(stats.me, align, [&](const Alignment &)
                {
                    const bool succeed = find(s.r_exons.begin(), s.r_exons.end(), align, f);
                    return options.filters.count(f.iID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                stats.ec.at(f.l)++;
                stats.ce.at(s.r_iso2Gene.at(f.iID))++;
            }
        }

        /*
         * Classify at the intron level
         */
        
        else
        {
            q_introns.push_back(align);

            if (classify(stats.mi, align, [&](const Alignment &)
                {
                    const bool succeed = checkSplice(align, f);
                    return options.filters.count(f.iID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                stats.ic.at(align.l)++;
                stats.ci.at(s.r_iso2Gene.at(f.iID))++;
            }
        }
    });
    
    /*
     * Classify at the base level
     */

    countBase(s.r_l_exons, q_exons, stats.mb, stats.cb);

    /*
     * Calculating for the references. The idea is similar to cuffcompare, where a reference is counted
     * for each true-positive. Anything that has not been detected will be the false-positives.
     */

    count_ref(stats.ec, stats.me.nr);
    count_ref(stats.ic, stats.mi.nr);
    stats.mb.nr = s.r_c_exons;

    assert(stats.me.nr && stats.mi.nr && stats.mb.nr);
    
    // The structure depends on the mixture
    const auto seqs = s.r_pair(options.mix);

    /*
     * Calculate for the LOS
     */

    stats.se = Expression::analyze(stats.ce, seqs);
    stats.si = Expression::analyze(stats.ci, seqs);
    stats.sb = Expression::analyze(stats.cb, seqs);

    AnalyzeReporter::report("ralign_base.stats", stats.mb, stats.sb, stats.cb, options.writer);
    AnalyzeReporter::report("ralign_exon.stats", stats.me, stats.se, stats.ce, options.writer);
    AnalyzeReporter::report("ralign_introns.stats", stats.mi, stats.si, stats.ci, options.writer);

	return stats;
}