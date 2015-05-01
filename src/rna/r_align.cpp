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

    auto cb = RAnalyzer::counter(Gene, options.mix);
    auto ce = RAnalyzer::counter(Gene, options.mix);
    auto cj = RAnalyzer::counter(Gene, options.mix);

    std::vector<Alignment> q_exons;
    std::vector<Alignment> q_juns;

    ParserSAM::parse(file, [&](const Alignment &align)
    {
        if (!align.mapped)
        {
            return;
        }

        Feature f;

        /*
         * Classify at the exon level
         */

        if (!align.spliced)
        {
            q_exons.push_back(align);

            if (classify(stats.me, align, [&](const Alignment &)
                {
                    return find(s.r_fs.begin(), s.r_fs.end(), align, f);
                }))
            {
                ce.at(s.r_iso2Gene.at(f.iID))++;
            }
        }

        /*
         * Classify at the junction level
         */
        
        else
        {
            Feature f1, f2;
            q_juns.push_back(align);
            
            if (classify(stats.mj, align, [&](const Alignment &)
                {
                    return checkSplice(s, align, f1, f2);
                }))
            {
                
                cj.at(s.r_iso2Gene.at(f1.iID))++;
            }
        }
    });
    
    /*
     * Classify at the base level
     */

    countBase(s.r_l_introns, q_juns,  stats.mb);
    countBase(s.r_l_exons,   q_exons, stats.mb);

    stats.me.nr() = q_exons.size();
    stats.mj.nr() = s.r_introns.size();
    stats.mb.nr() = s.r_c_exons + s.r_c_introns;

    // The structure depends on the mixture
    const auto seqs = s.r_pair(options.mix);

    /*
     * Calculate for the sensitivity
     */

    stats.se = Expression::analyze(ce, seqs);
    stats.sj = Expression::analyze(cj, seqs);
    stats.sb = Expression::analyze(cb, seqs);

    const auto &writer = options.writer;
    
    // Report for the base-level
    AnalyzeReporter::reportClassify("ralign_base.stats", stats.mb, stats.sb, cb, writer);
    
    // Report for the exon-level
    AnalyzeReporter::reportClassify("ralign_exon.stats", stats.me, stats.se, ce, writer);

    // Report for the intron-level
    AnalyzeReporter::reportClassify("ralign_junction.stats", stats.mj, stats.sj, cj, writer);

	return stats;
}