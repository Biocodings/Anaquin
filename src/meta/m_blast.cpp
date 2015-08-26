#include <numeric>
#include "meta/m_blast.hpp"
#include <boost/format.hpp>
#include "parsers/parser_blast.hpp"

using namespace Anaquin;

MBlast::Stats MBlast::stats(const FileName &file, const Options &options)
{
    /*
     * Create data-strucutre for the sequins
     */
    
    SequinAlign m;

    for (const auto &id : Standard::instance().r.seqs())
    {
        const auto &d = Standard::instance().r.seq(id, MixA);

        m[id].id   = id;
        m[id].seqA = &d;

        //if (mixB.count(m[seqID].id))
        //{
        //    m[seqID].seqB = mixB.at(m[seqID].id);
        //}
    }
/*
    const auto &mixB = Standard::instance().seqs_2;
    
    for (const auto &seq : Standard::instance().seqs_1)
    {
        const auto &seqID = seq.first;
        
        m[seqID].id   = seqID;
        m[seqID].seqA = seq.second;
        
        if (mixB.count(m[seqID].id))
        {
            m[seqID].seqB = mixB.at(m[seqID].id);
        }
    }
*/
    /*
     * Create data-strucutre for the alignment
     */
    
    options.info("Comparing alignment");
    
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
            m.at(id).contigs.push_back(contig);
        }
        else
        {
            options.warn((boost::format("%1% is not a sequin") % id).str());
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
            
            Base gaps = 0;
            Base match = 0;
            Base mismatch = 0;
            
            for (const auto &contig : align.contigs)
            {
                gaps     += contig.gap;
                match    += contig.match;
                mismatch += contig.mismatch;
            }
            
            //assert(align.seqA.length == align.seqB.length);
            assert(match > mismatch && match > gaps);
            
            // Fraction of sequins covered by alignments
            align.covered = (double) total / align.seqA->length;
            
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

void MBlast::analyze(const std::string &file, const Options &options)
{
    const auto stats = MBlast::stats(file);

    options.info("Generating statistics");

    /*
     * Generate summary statistics
     */

    {
        options.writer->open("MetaPSL_summary.stats");
        
        const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        options.writer->write((boost::format(format) % "id"
                                                     % "contigs"
                                                     % "cover"
                                                     % "mismatch"
                                                     % "gap").str());
        
        for (const auto &i : stats.metas)
        {
            const auto &align = i.second;
            
            options.writer->write((boost::format(format) % align.id
                                                         % align.contigs.size()
                                                         % align.covered
                                                         % align.mismatch
                                                         % align.gaps).str());
        }

        options.writer->close();
    }
}