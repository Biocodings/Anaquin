#include <numeric>
#include "meta/m_blat.hpp"
#include <boost/format.hpp>
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

MBlat::Stats MBlat::analyze(const FileName &file, const Options &o)
{
    /*
     * Create data-structure for each of the sequin
     */
    
    SequinAlign m;

    for (const auto &i : Standard::instance().r_meta.data())
    {
        m[i.second.id] = std::shared_ptr<MetaAlignment>(new MetaAlignment());
        m[i.second.id]->seq = &i.second;
    }

    MBlat::Stats stats;

    /*
     * Create data-structure for the alignment
     */
    
    ParserBlast::parse(file, [&](const ParserBlast::BlastLine &l, const ParserProgress &)
    {
        // Eg: M2_G, M10_G
        const SequinID id = l.tName;

        if (m.count(id))
        {
            stats.n_chrT++;

            AlignedContig contig;
            
            contig.id       = l.qName;
            contig.l        = Locus(l.tStart, l.tEnd);
            contig.match    = l.matches;
            contig.mismatch = l.mismatch;
            
            // Only interested in the target (eg: M10_G)
            contig.gap = l.tGaps;
            
            // That's because we might have multiple contigs aligned to a sequin
            m.at(id)->contigs.push_back(contig);
        }
        else
        {
            stats.n_hg38++;
            o.warn((boost::format("%1% is not a sequin") % id).str());
        }
    });

    /*
     * Traverse through all sequins, and calculate statistics for all alignments for each of those sequin.
     */
    
    // For each sequin in the reference...
    for (auto &i : m)
    {
        const auto &id = i.first;
        
        // Aligments for this particular sequin
        auto &align = i.second;

        if (!align->contigs.empty())
        {
            // Make sure the contigs are non-overlapping
            const auto merged = Locus::merge<AlignedContig, Locus>(align->contigs);

            // Total non-overlapping bases for the alignments
            const auto total = std::accumulate(merged.begin(), merged.end(), 0, [&](int sum, const Locus &l)
            {
                return sum + l.length();
            });

            /*
             * Don't consider for overlapping because a base can be matched or mismatched.
             */
            
            Base gaps     = 0;
            Base sums     = 0;
            Base match    = 0;
            Base mismatch = 0;
            
            /*
             * Generating summary statistic for this sequin
             */
            
            for (const auto &contig : align->contigs)
            {
                gaps     += contig.gap;
                match    += contig.match;
                mismatch += contig.mismatch;
                sums     += (gaps + match + mismatch);
            }

            assert(align->seq->length);

            // Fraction of sequins covered by alignment
            align->covered = static_cast<double>(total) / align->seq->length;
            
            // Fraction of mismatch bases in alignment
            align->mismatch = static_cast<double>(mismatch) / match;
            
            // Fraction of bases covered by gaps
            align->gaps = static_cast<double>(gaps) / sums;

            stats.oGaps     += gaps;
            stats.oMatch    += match;
            stats.oMismatch += mismatch;
            stats.total     += align->seq->length;
            
            if (align->covered  > 1) { o.warn((boost::format("%1% coverage: %2%") % id % align->covered).str()); }
            if (align->mismatch > 1) { o.warn((boost::format("%1% mismatch: %2%") % id % align->mismatch).str()); }

            assert(align->gaps    >= 0.0 && align->gaps     <= 1.0);
            assert(align->covered >= 0.0 && align->mismatch >= 0.0);
            
            // Create an alignment for each contig that aligns to the MetaQuin
            for (const auto &i : align->contigs)
            {
                stats.aligns[i.id] = align;
            }
        }

        stats.metas[align->seq->id] = align;
    }

    return stats;
}

void MBlat::report(const FileName &file, const Options &o)
{
    const auto stats = MBlat::analyze(file);

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    {
        o.writer->open("MetaPSL_summary.stats");
        
        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Community: %2%\n"
                             "   Synthetic: %3%\n\n"
                             "   Contigs: %4%\n"
                             "   Assembled: %5%\n"
                             "   Reference: %6%\n\n"
                             "   ***\n"
                             "   *** The following overlapping statistics are computed by proportion\n"
                             "   ***\n\n"
                             "   Match: %7%\n"
                             "   Gaps: %8%\n"
                             "   Mismatch: %9%\n";

        o.writer->write((boost::format(summary) % file
                                                % stats.n_hg38
                                                % stats.n_chrT
                                                % stats.aligns.size()
                                                % stats.sequin()
                                                % stats.metas.size()
                                                % stats.overMatch()
                                                % stats.overGaps()
                                                % stats.overMismatch()).str());
        o.writer->close();
    }

    {
        o.writer->open("MetaPSL_quins.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.writer->write((boost::format(format) % "id"
                                               % "contigs"
                                               % "covered"
                                               % "mismatch"
                                               % "gaps").str());
        
        for (const auto &i : stats.metas)
        {
            const auto &align = i.second;
            
            o.writer->write((boost::format(format) % align->seq->id
                                                   % align->contigs.size()
                                                   % align->covered
                                                   % align->mismatch
                                                   % align->gaps).str());
        }

        o.writer->close();
    }
}