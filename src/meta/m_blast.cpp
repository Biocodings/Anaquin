#include <numeric>
#include <boost/format.hpp>
#include "meta/m_blast.hpp"
#include "data/standard.hpp"
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

MBlast::Stats MBlast::analyze(const std::string &file, const Options &options)
{
    /*
     * Create data-strucutre for each sequin
     */
    
    options.info("Loading mixture file");
    
    std::map<SequinID, MetaAlignment> m;

    const auto &mixB = Standard::instance().seqs_2;

    for (const auto &seq : Standard::instance().seqs_1)
    {
        const auto &seqID = seq.first;

        m[seqID].id   = seq.first;
        m[seqID].seqA = seq.second;
        
        if (mixB.count(m[seqID].id))
        {
            m[seqID].seqB = mixB.at(m[seqID].id);
        }
    }

    options.info("Comparing alignment with sequins");

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

            m.at(id).contigs.push_back(contig);
        }
        else
        {
            options.warn((boost::format("%1% is not a sequin (given in alignment)") % id).str());
        }
    });

    /*
     * The variable m now holds information on the aligned contigs.
     * Traverse through all sequins, and calculate statistics for all alignments
     * to each of those sequin.
     */

    MBlast::Stats stats;

    for (auto &i : m)
    {
        // The data-structure for the aligment for this particular sequin
        auto &align = i.second;

        if (!align.contigs.empty())
        {
            // Make sure the contigs are non-overlapping
            const auto merged = Locus::merge<AlignedContig, Locus>(align.contigs);

            // The total non-overlapping bases for the alignments
            const auto total = std::accumulate(merged.begin(), merged.end(), 0, [&](int sum, const Locus &l)
            {
                return sum + l.length();
            });

            /*
             * Don't consider for overlapping because a base can be matched or mismatched.
             */
            
            BasePair gaps = 0;
            BasePair match = 0;
            BasePair mismatch = 0;

            for (const auto &contig : align.contigs)
            {
                gaps     += contig.gap;
                match    += contig.match;
                mismatch += contig.mismatch;
            }
            
            //assert(align.seqA.length == align.seqB.length);
            assert(match > mismatch && match > gaps);
            
            // Fraction of sequins covered by alignments
            align.covered = (double) total / align.seqA.length;
            
            // Fraction of mismatch bases in alignments
            align.mismatch = (double) mismatch / match;

            // Fraction of bases covered by gaps
            align.gaps = (double) gaps / match;

            assert(align.gaps     >= 0.0 && align.gaps     <= 1.0);
            assert(align.covered  >= 0.0 && align.covered  <= 1.0);
            assert(align.mismatch >= 0.0 && align.mismatch <= 1.0);

            // Create an alignment for each contig that aligns to the metaquin
            for (const auto &i : align.contigs)
            {
                stats.aligns[i.id] = align;
            }
        }

        stats.metas[align.id] = align;
    }

    options.info("Generating statistics");

    {
        options.writer->open("MetaPSL_summary.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        options.writer->write((boost::format(format) % "id"
                                                     % "contigs"
                                                     % "cover"
                                                     % "mismatch"
                                                     % "gap").str());
        
        for (const auto &align : stats.metas)
        {
            options.writer->write((boost::format(format) % align.second.id
                                                         % align.second.contigs.size()
                                                         % align.second.covered
                                                         % align.second.mismatch
                                                         % align.second.gaps).str());
        }
        
        options.writer->close();
    }

    return stats;
}