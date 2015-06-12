#include <map>
#include <numeric>
#include <iostream>
#include "meta/m_blast.hpp"
#include "data/standard.hpp"
#include "parsers/parser_blast.hpp"

using namespace Spike;

MBlast::Stats MBlast::analyze(const std::string &file)
{
    std::map<std::string, ParserBlast::BlastLine> psl;

    ParserBlast::parse(file, [&](const ParserBlast::BlastLine &l, const ParserProgress &)
    {
        psl[l.qName] = l;
    });

    if (psl.empty())
    {
        throw EmptyFileError(file);
    }

    /*
     * Create data-strucutre for each metaquin
     */
    
    std::map<MetaQuinID, MetaAlignment> m;

    const auto &mixB = Standard::instance().m_seq_B;
    
    for (const auto &seq : Standard::instance().m_seq_A)
    {
        m[seq.first].id   = seq.first;
        m[seq.first].seqA = seq.second;
        m[seq.first].seqB = mixB.at(m[seq.first].id);
    }
    
    /*
     * Compare each alignment to the metaquins
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

            m.at(id).contigs.push_back(contig);
        }
        else
        {
            std::cout << "Warning: " << id << " not a metaquin (given in alignment)" << std::endl;
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
            
            assert(match > mismatch && match > gaps);
            assert(align.seqA.l.length() == align.seqB.l.length());
            
            // Fraction of sequins covered by alignments
            align.covered = (double) total / align.seqA.l.length();
            
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

    return stats;
}