#include <numeric>
#include "meta/m_blat.hpp"
#include <boost/format.hpp>
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

MBlast::Stats MBlast::stats(const FileName &file, const Options &o)
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

    /*
     * Create data-structure for the alignment
     */
    
    ParserBlast::parse(file, [&](const ParserBlast::BlastLine &l, const ParserProgress &)
    {
        // Eg: M2_G, M10_G
        const SequinID id = l.tName;

        if (m.count(id))
        {
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
            o.warn((boost::format("%1% is not a sequin") % id).str());
        }
    });

    /*
     * Traverse through all sequins, and calculate statistics for all alignments
     * for each of those sequin.
     */
    
    MBlast::Stats stats;

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
            
            Base gaps  = 0;
            Base sums  = 0;
            Base match = 0;
            Base mismatch = 0;
            
            for (const auto &contig : align->contigs)
            {
                gaps     += contig.gap;
                match    += contig.match;
                mismatch += contig.mismatch;
                sums     += (gaps + match + mismatch);
            }

            // Fraction of sequins covered by alignments
            align->covered = (double) total / align->seq->length;

            // Fraction of mismatch bases in alignments
            align->mismatch = (double) mismatch / match;
            
            // Fraction of bases covered by gaps
            align->gaps = (double) gaps / sums;

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

void MBlast::report(const std::string &file, const Options &o)
{
    const auto stats = MBlast::stats(file);

    o.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    {
        o.writer->open("MetaPSL_summary.stats");
        
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