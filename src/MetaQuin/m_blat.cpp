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
    const auto r1 = r.regs1();
    
    for (const auto &i : r.seqsL1())
    {
        m[i] = std::shared_ptr<MetaAlignment>(new MetaAlignment());
        m[i]->seq.id = i;
        m[i]->seq.l  = Locus(1, r1.at(i).length());
    }
    
    MBlat::Stats stats;
    
    /*
     * Create data-structure for the alignments
     */
    
    ParserBlat::parse(file, [&](const ParserBlat::Data &x)
    {
        // Eg: MQ_16
        const auto id = x.tName;
        
        if (m.count(id))
        {
            stats.nSeqs++;
            
            AlignedContig c;
            
            c.id = x.qName;
            
            // (c.l != c.qSize) because a contig can be mapped to a sequin multiple times
            c.l  = Locus(x.tStart, x.tEnd);

            c.match    = x.match;
            c.mismatch = x.mismatch;
            
            c.rGap   = x.tGap;
            c.rStart = x.tStart;
            c.rEnd   = x.tEnd;
            c.rSize  = x.tSize;
            
            c.qGap   = x.qGap;
            c.qStart = x.qStart;
            c.qEnd   = x.qEnd;
            c.qSize  = x.qSize;
            
            c.qGapCount = x.qGapCount;
            c.rGapCount = x.tGapCount;
            
            // That's because we might have multiple contigs aligned to a sequin
            m.at(id)->contigs.push_back(c);
            
            A_ASSERT(c.l.length());
            
            stats.t2l[x.tName] = x.tSize;
            
            /*
             * Building mappings for contigs
             */
            
            // Size of the contig
            stats.c2l[c.id] = c.qSize;
            
            // Size of the aligned contig (<= contig size)
            stats.c2a[c.id] = c.qSize; // c.l.length() - 1;
            
            A_ASSERT(stats.c2a[c.id] <= stats.c2l[c.id]);
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
            
            if (align->seq.l.length() > 2)
            {
                l = align->seq.l.length();
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
            stats.total     += align->seq.l.length();
            
            if (align->oRGaps > 1)    { o.warn((boost::format("%1% (ga): %2%") % id % align->oRGaps).str());    }
            if (align->oQGaps > 1)    { o.warn((boost::format("%1% (ga): %2%") % id % align->oQGaps).str());    }
            if (align->covered > 1)   { o.warn((boost::format("%1% (co): %2%") % id % align->covered).str());   }
            if (align->oMismatch > 1) { o.warn((boost::format("%1% (mm): %2%") % id % align->oMismatch).str()); }
            
            //assert(align->oGaps   >= 0.0 && align->oGaps     <= 1.0);
            assert(align->covered >= 0.0 && align->oMismatch >= 0.0);
            
            // Create an alignment for each contig that aligns to the MetaQuin
            for (const auto &i : align->contigs)
            {
                stats.c2s[i.id] = align->seq.id;
                stats.aligns[i.id] = align;
            }
        }
        
        stats.metas[align->seq.id] = align;
    }
    
    return stats;
}
