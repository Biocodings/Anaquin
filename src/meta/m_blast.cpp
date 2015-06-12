#include <map>
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

            contig.id = l.qName;
            contig.l  = Locus(l.tStart, l.tEnd);

            m.at(id).contigs.push_back(contig);
        }
        else
        {
            std::cout << "Warning: " << id << " not a metaquin" << std::endl;
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
            
            BasePair gaps  = 0;
            BasePair total = 0;
            
            for (std::size_t j = 0; j < merged.size(); j++)
            {
                if (j && merged[j-1].end + 1 != merged[j].start)
                {
                    gaps += Locus(merged[j-1].end + 1, merged[j].start - 1).length();
                }
                
                total += merged[j].length();
            }
            
            assert(align.seqA.l.length() == align.seqB.l.length());
            
            // Fraction of bases covered by alignments
            align.covered = (double) total / align.seqA.l.length();
            
            // Fraction of bases not covered by alignments
            align.mismatch = 1 - align.covered;
            
            // Fraction of bases covered by gaps
            align.gaps = (double) gaps / align.seqA.l.length();
            
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