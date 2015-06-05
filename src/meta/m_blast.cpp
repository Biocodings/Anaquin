#include <map>
#include <iostream>
#include "standard.hpp"
#include "meta/m_blast.hpp"
#include "parsers/parser_blat.hpp"

using namespace Spike;

MBlast::Stats MBlast::analyze(const std::string &file)
{
    std::map<std::string, ParserBlat::BlatLine> psl;

    ParserBlat::parse(file, [&](const ParserBlat::BlatLine &l, const ParserProgress &)
    {
        psl[l.qName] = l;
    });

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
    
    ParserBlat::parse(file, [&](const ParserBlat::BlatLine &l, const ParserProgress &)
    {
        if (m.count(l.tName))
        {
            m.at(l.tName).ids.push_back(l.tName);
            m.at(l.tName).aligns.push_back(Locus(l.tStart, l.tEnd));
        }
        else
        {
            std::cout << "Warning: " << l.tName << " not a metaquin" << std::endl;
        }
    });
    
    /*
     * Convert the results
     */
    
    MBlast::Stats stats;

    for (auto &i : m)
    {
        /*
         * In order to calculate for the coverage, we'll need to merge the contigs. Although they
         * are likely already non-overlapping, just to make sure...
         */
        
        const auto merged = Locus::merge<Locus, Locus>(i.second.aligns);

        BasePair gaps  = 0;
        BasePair total = 0;
        
        for (std::size_t j = 0; j < merged.size(); j++)
        {
            // Only calculate if there is indeed a gap (although almost always there is)
            if (j && merged[j-1].end + 1 != merged[j].start)
            {
                gaps += Locus(merged[j-1].end + 1, merged[j].start - 1).length();
            }
            
            total += merged[j].length();
        }
        
        assert(i.second.seqA.l.length() == i.second.seqB.l.length());
        
        // Fraction of bases covered by alignments
        i.second.coverage = (double) total / i.second.seqA.l.length();
        
        // Fraction of bases not covered by alignments
        i.second.mismatch = 1 - i.second.coverage;
        
        // Fraction of bases covered by gaps
        i.second.gaps = (double) gaps / i.second.seqA.l.length();
        
        assert(i.second.gaps     >= 0.0 && i.second.gaps     <= 1.0);
        assert(i.second.coverage >= 0.0 && i.second.coverage <= 1.0);
        assert(i.second.mismatch >= 0.0 && i.second.mismatch <= 1.0);

        // Create an alignment for each contig that aligns to the metaquin
        for (const ContigID &id : i.second.ids)
        {
            stats.aligns[id] = i.second;
        }

        stats.metas.insert(i.second);
    }

    return stats;
}