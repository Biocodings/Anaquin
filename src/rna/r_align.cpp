#include "rna/r_align.hpp"
#include "stats/expression.hpp"
#include "parsers/parser_sam.hpp"

using namespace Spike;

// Find the matching intron by locus given a spliced alignment
static bool findIntron(const Alignment &align, Feature &f)
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

void RAlign::reportGeneral(const std::string &file, const RAlign::Stats &stats, const Options &options)
{
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

    options.writer->open(file);
    options.writer->write((boost::format(format) % "samples"
                                                 % "chrT"
                                                 % "dilution"
                                                 % "sn"
                                                 % "sp"
                                                 % "spliced sn"
                                                 % "spliced sp"
                                                 % "los").str());
    options.writer->write((boost::format(format) % stats.n_samps
                                                 % stats.n_chrT
                                                 % stats.dilution()
                                                 //% stats.
/*
                                                 % ss.n_seqs
                                                 % ss.n_samps
                                                 % ss.dilution()
 */
 
 ).str());
    options.writer->close();
}

RAlign::Stats RAlign::analyze(const std::string &file, const Options &options)
{
    RAlign::Stats stats;
    const auto &s = Standard::instance();

    std::vector<Alignment> exons, introns;

    options.terminal->write("Parsing alignments");
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if ((p.i % 1000) == 0)
        {
            options.terminal->write("Processed: " + std::to_string(p.i));
        }

        Feature f;

        if (!align.mapped)
        {
            return;
        }
        else if (align.id != s.id)
        {
            stats.n_samps++;
            return;
        }

        stats.n_chrT++;
        
        // Whether the read has mapped to a feature correctly
        bool succeed = false;

        options.logger->write((boost::format("%1% %2% %3%") % align.id % align.l.start % align.l.end).str());

        /*
         * Collect statistics at the exon level
         */

        if (!align.spliced)
        {
            exons.push_back(align);

            if (classify(stats.pe.m, align, [&](const Alignment &)
                {
                    succeed = find(s.r_exons.begin(), s.r_exons.end(), align, f);
                    return options.filters.count(f.tID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                stats.ce.at(s.r_isoformToGene.at(f.tID))++;
            }
        }

        /*
         * Collect statistics at the intron level
         */

        else
        {
            introns.push_back(align);

            if (classify(stats.pi.m, align, [&](const Alignment &)
                {
                    succeed = findIntron(align, f);
                    return options.filters.count(f.tID) ? Ignore : succeed ? Positive : Negative;
                }))
            {
                stats.ci.at(s.r_isoformToGene.at(f.tID))++;
            }
        }

        /*
         * Collect statistics at the gene level
         */
        
        const auto geneID = f.tID;

        if (succeed && stats.g_exon_tracker.count(geneID))
        {
            if (align.spliced)
            {
                auto &p = stats.g_intron_tracker.at(geneID);

                classify(p.m, align, [&](const Alignment &)
                {
                    // This is only executed if succeed is true
                    return Positive;
                });
            }
            else
            {
                auto &p = stats.g_exon_tracker.at(geneID);

                classify(p.m, align, [&](const Alignment &)
                {
                    // This is only executed if succeed is true
                    return Positive;
                });
            }
        }
    });

    options.terminal->write("Counting references");
    
    /*
     * Calculate for references. The idea is similar to cuffcompare, each true-positive is counted
     * as a reference. Anything that is undetected in the experiment will be counted as a single reference.
     */

    sums(stats.ce, stats.pe.m.nr);
    sums(stats.ci, stats.pi.m.nr);

    options.terminal->write("Counting bases");
    
    /*
     * The counts for query bases is the total non-overlapping length of all the exons in the experiment.
     * The number is expected to approach the reference length (calculated next) for a very large
     * experiment with sufficient coverage.
     */

    countBase(s.r_l_exons, exons, stats.pb.m, stats.cb);

    /*
     * The counts for references is the total length of all known non-overlapping exons.
     * For example, if we have the following exons:
     *
     *    {1,10}, {50,55}, {70,74}
     *
     * The length of all the bases is 10+5+4 = 19.
     */
    
    stats.pb.m.nr = s.r_c_exons;

    assert(stats.pe.m.nr && stats.pi.m.nr && stats.pb.m.nr);

    // The structure depends on the mixture
    const auto seqs = s.r_gene(options.mix);

    /*
     * Calculate for the LOS
     */

    stats.pe.s = Expression::analyze(stats.ce, seqs);
    stats.pi.s = Expression::analyze(stats.ci, seqs);
    stats.pb.s = Expression::analyze(stats.cb, seqs);

    options.terminal->write("Writing results");
    
    // Write out general statistics
    //reportGeneral("ralign_general.statsD", stats, options);

    /*
     * Write out statistics for various levels
     */

    AnalyzeReporter::report("ralign_base.stats",    stats.pb, stats.cb, options.writer);
    AnalyzeReporter::report("ralign_exon.stats",    stats.pe, stats.ce, options.writer);
    AnalyzeReporter::report("ralign_introns.stats", stats.pi, stats.ci, options.writer);

	return stats;
}