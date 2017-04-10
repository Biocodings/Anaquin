#include <numeric>
#include <boost/format.hpp>
#include "MetaQuin/m_blat.hpp"
#include "parsers/parser_blat.hpp"

using namespace Anaquin;

MBlat::Stats MBlat::analyze(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_meta;
    
    /*
     * Create data-structure for each of the sequin
     */
    
    SequinAlign m;
    
    for (const auto &i : r.data())
    {
        m[i.second.id] = std::shared_ptr<MetaAlignment>(new MetaAlignment());
        m[i.second.id]->seq = &i.second;
    }
    
    MBlat::Stats stats;
    
    /*
     * Create data-structure for the alignments
     */
    
    ParserBlat::parse(file, [&](const ParserBlat::Data &l, const ParserProgress &)
    {
        // Eg: M2_G, M10_G
        const auto id = l.tName;
        
        if (m.count(id))
        {
            stats.nSeqs++;
            
            AlignedContig contig;
            
            contig.id = l.qName;
            contig.l  = Locus(l.tStart, l.tEnd);
            
            contig.match    = l.match;
            contig.mismatch = l.mismatch;
            
            contig.rGap   = l.tGap;
            contig.rStart = l.tStart;
            contig.rEnd   = l.tEnd;
            contig.rSize  = l.tSize;
            
            contig.qGap   = l.qGap;
            contig.qStart = l.qStart;
            contig.qEnd   = l.qEnd;
            contig.qSize  = l.qSize;
            
            contig.qGapCount = l.qGapCount;
            contig.rGapCount = l.tGapCount;
            
            // That's because we might have multiple contigs aligned to a sequin
            m.at(id)->contigs.push_back(contig);
            
            A_ASSERT(contig.l.length());
            
            /*
             * TODO: What about a contig being mapped to the same sequin multiple times?
             */
            
            stats.t2l[l.tName] = l.tSize;
            
            /*
             * Building mappings for contigs
             */
            
            // Size of the contig
            stats.c2l[contig.id] = contig.qSize;
            
            // Size of the aligned contig (<= contig size)
            stats.c2a[contig.id] = contig.l.length() - 1;
            
            A_ASSERT(stats.c2a[contig.id] <= stats.c2l[contig.id]);
        }
        else
        {
            stats.nEndo++;
            o.warn((boost::format("%1% is not a sequin") % id).str());
        }
    });
    
    /*
     * Traverse through the sequins, and calculate statistics for all alignments for each of those sequin.
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
            
            /*
             * Generating non-overlapping summary statistic for this sequin
             */
            
            const auto total = std::accumulate(merged.begin(), merged.end(), 0, [&](int sum, const Locus &l)
            {
                return sum + l.length();
            });
            
            /*
             * Generating overlapping statistics for this sequin. In this context, target refers
             * to the sequin whereas query refers to the contig.
             */
            
            Base oRGaps    = 0;
            Base oQGaps    = 0;
            Base oMatch    = 0;
            Base oMismatch = 0;
            Base qSums     = 0;
            
            for (const auto &contig : align->contigs)
            {
                qSums     += contig.qSize;
                oRGaps    += contig.rGap;
                oQGaps    += contig.qGap;
                oMatch    += contig.match;
                oMismatch += contig.mismatch;
            }
            
            /*
             * We need the sequin length. Try to get it from the BED annotation. Otherwise,
             * try from the alignment file.
             */
            
            // Sequin length
            Base l;
            
            if (align->seq->l.length() > 2)
            {
                l = align->seq->l.length();
            }
            else if (stats.t2l.count(align->id()))
            {
                l = stats.t2l[align->id()];
            }
            else
            {
                throw std::runtime_error("Sequin length for " + align->id() + " not found");
            }
            
            A_ASSERT(l > 2);
            
            // Proportion of non-overlapping bases covered or assembled
            align->covered = static_cast<double>(total) / l;
            
            // Proportion of overlapping matches relative to the sequin
            align->oMatch = static_cast<double>(oMatch) / l;
            
            // Proportion of overlapping mismatches relative to the sequin
            align->oMismatch = static_cast<double>(oMismatch) / l;
            
            // Proportion of overlapping gaps relative to the sequin
            align->oRGaps = static_cast<double>(oRGaps) / l;
            
            // Proportion of overlapping gaps relative to the sequin
            align->oQGaps = static_cast<double>(oQGaps) / qSums;
            
            /*
             * Update overall statistics for the sequins
             */
            
            // Overall total gaps in all sequins
            stats.oRGaps += oRGaps;
            
            // Overall total gaps in all contigs
            stats.oQGaps += oQGaps;
            
            // Overall total size in all contigs
            stats.oQSums += qSums;
            
            stats.oMatch += oMatch;
            stats.oMismatch += oMismatch;
            stats.total     += align->seq->l.length();
            
            if (align->oRGaps > 1)    { o.warn((boost::format("%1% (ga): %2%") % id % align->oRGaps).str());    }
            if (align->oQGaps > 1)    { o.warn((boost::format("%1% (ga): %2%") % id % align->oQGaps).str());    }
            if (align->covered > 1)   { o.warn((boost::format("%1% (co): %2%") % id % align->covered).str());   }
            if (align->oMismatch > 1) { o.warn((boost::format("%1% (mm): %2%") % id % align->oMismatch).str()); }
            
            //assert(align->oGaps   >= 0.0 && align->oGaps     <= 1.0);
            assert(align->covered >= 0.0 && align->oMismatch >= 0.0);
            
            // Create an alignment for each contig that aligns to the MetaQuin
            for (const auto &i : align->contigs)
            {
                stats.c2s[i.id] = align->seq->id;
                stats.aligns[i.id] = align;
            }
        }
        
        stats.metas[align->seq->id] = align;
    }
    
    return stats;
}
