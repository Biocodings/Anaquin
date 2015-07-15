#include <numeric>
#include <iostream>
#include <boost/format.hpp>
#include "meta/m_blast.hpp"
#include "data/standard.hpp"
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

MBlast::Stats MBlast::analyze(const std::string &file, const AnalyzerOptions &options)
{
    std::map<std::string, ParserBlast::BlastLine> psl;

    ParserBlast::parse(file, [&](const ParserBlast::BlastLine &l, const ParserProgress &)
    {
        psl[l.qName] = l;
    });

    if (psl.empty())
    {
        throw InvalidFileError(file);
    }

    /*
     * Create data-strucutre for each metaquin
     */
    
    std::map<MetaQuinID, MetaAlignment> m;

    const auto &mixB = Standard::instance().m_seqs_B;
    
    for (const auto &seq : Standard::instance().m_seqs_A)
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
            
            assert(align.seqA.length == align.seqB.length);
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

    /*
     * Write out results
     */

    options.writer->open("align_stats.txt");
    
    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";

    options.writer->write((boost::format(format) % "ID"
                                                 % "ConA"
                                                 % "ConB"
                                                 % "Contigs"
                                                 % "Covered"
                                                 % "Mismatch"
                                                 % "Gaps").str());

    for (const auto &align : stats.metas)
    {
        options.writer->write((boost::format(format)
                                    % align.second.id
                                    % align.second.seqA.abund()
                                    % align.second.seqB.abund()
                                    % align.second.contigs.size()
                                    % align.second.covered
                                    % align.second.mismatch
                                    % align.second.gaps).str());
        options.writer->write("\n");
    }

    options.writer->close();

    return stats;
}