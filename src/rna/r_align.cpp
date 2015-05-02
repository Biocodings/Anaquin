#include <iostream>
#include <assert.h>
#include "r_align.hpp"
#include "biology.hpp"
#include "expression.hpp"
#include <boost/format.hpp>
#include "writers/writer.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

static bool checkSplice(const Standard &s, const Alignment &align, Feature &e1, Feature &e2)
{
    assert(align.spliced);

    /*
     * A spliced alignment is correct if it aligns to two consecutive exons
     */
        
    for (auto i = 0; i < s.r_exons.size(); i++)
    {
        e1 = s.r_exons[i];

        // Check if it starts inside exon_1
        if (align.l.start >= e1.l.start && align.l.start <= e1.l.end && align.l.end > e1.l.end
            && i != s.r_exons.size() - 1)
        {
            for (auto j = i + 1; j < s.r_exons.size(); j++)
            {
                e2 = s.r_exons[j];

                if (e1.iID == e2.iID)
                {
                    // Check if it ends inside exon_2
                    if (align.l.start < e2.l.start && align.l.end >= e2.l.start && align.l.end <= e2.l.end)
                    {
                        return true;
                    }
                }
            }
        }
    }
    
    return false;
}

RAlignStats RAlign::analyze(const std::string &file, const Options &options)
{
    RAlignStats stats;
    const auto &s = Standard::instance();

    stats.cb = RAnalyzer::counter(Gene, options.mix);
    stats.ce = RAnalyzer::counter(Gene, options.mix);
    stats.ci = RAnalyzer::counter(Gene, options.mix);

    std::vector<Alignment> q_exons, q_introns;

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &)
    {
        if (!align.mapped)
        {
            return;
        }

        /*
         * Classify at the exon level
         */

        if (!align.spliced)
        {
            Feature f;
            q_exons.push_back(align);

            if (classify(stats.me, align, [&](const Alignment &)
                {
                    const bool succeed = find(s.r_exons.begin(), s.r_exons.end(), align, f);
                    return options.filters.count(f.iID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                stats.ce.at(s.r_iso2Gene.at(f.iID))++;
            }
        }

        /*
         * Classify at the intron level
         */
        
        else
        {
            Feature f1, f2;
            q_introns.push_back(align);

            if (classify(stats.mi, align, [&](const Alignment &)
                {
                    const bool succeed = checkSplice(s, align, f1, f2);
                    return options.filters.count(f1.iID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                assert(f1.iID == f2.iID);
                stats.ci.at(s.r_iso2Gene.at(f1.iID))++;
            }
        }
    });

    /*
     * Classify at the base level
     */

    countBase(s.r_l_exons,   q_exons, stats.mb);
    countBase(s.r_l_introns, q_introns, stats.mb);

    stats.me.nr = s.r_exons.size();
    stats.mi.nr = s.r_introns.size();
    stats.mb.nr = s.r_c_exons + s.r_c_introns;

    // The structure depends on the mixture
    const auto seqs = s.r_pair(options.mix);

    /*
     * Calculate for the sensitivity
     */

    stats.se = Expression::analyze(stats.ce, seqs);
    stats.si = Expression::analyze(stats.ci, seqs);
    stats.sb = Expression::analyze(stats.cb, seqs);

    const auto &writer = options.writer;
    
    // Report for the base-level
    AnalyzeReporter::report("ralign_base.stats", stats.mb, stats.sb, stats.cb, writer);
    
    // Report for the exon-level
    AnalyzeReporter::report("ralign_exon.stats", stats.me, stats.se, stats.ce, writer);

    // Report for the intron-level
    AnalyzeReporter::report("ralign_junction.stats", stats.mi, stats.si, stats.ci, writer);

	return stats;
}