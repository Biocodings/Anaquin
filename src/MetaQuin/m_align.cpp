#include "tools/tools.hpp"
#include "MetaQuin/m_align.hpp"
#include "MetaQuin/MetaQuin.hpp"
#include "parsers/parser_bam.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace Anaquin;

static MAlign::Stats init()
{
    const auto &r = Standard::instance().r_meta;
    
    MAlign::Stats stats;
    
    stats.inters = r.mInters();
    A_ASSERT(!stats.inters.empty());
    
    for (const auto &i : stats.inters)
    {
        const auto &cID = i.first;
        
        /*
         * We'd like to know the length of the genome but we don't have the information.
         * It doesn't matter because we can simply initalize it to the the maximum possible.
         */
        
        stats.data[cID].bLvl.fp = std::shared_ptr<MergedInterval>
                (new MergedInterval(cID, Locus(1, std::numeric_limits<Base>::max())));
    }
    
    return stats;
}

static void classifyAlign(MAlign::Stats &stats, ParserBAM::Data &align)
{
    if (!stats.data.count(align.cID))
    {
        return;
    }
    
    auto &x = stats.data.at(align.cID);
    
    Locus l;
    bool spliced;
    
    bool isTP = true;
    
    while (align.nextCigar(l, spliced))
    {
        Base lGaps = 0, rGaps = 0;
        bool isContained = false;
        
        auto f = [&](MergedInterval *m)
        {
            m->map(l, &lGaps, &rGaps);
            
            if (isContained)
            {
                lGaps = rGaps = 0;
            }
            
            const auto covered = (l.length() - lGaps - rGaps);
            
            stats.data[align.cID].lGaps[m->name()] += lGaps;
            stats.data[align.cID].lGaps[m->name()] += rGaps;
            stats.data[align.cID].align[m->name()] += covered;
            
            A_ASSERT(covered >= 0);
            A_ASSERT(l.length() > lGaps);
            A_ASSERT(l.length() > rGaps);
        };
        
        // Does the read aligned within a region (eg: gene)?
        const auto m = stats.inters.at(align.cID).contains(l);
        
        if (m)
        {
            isContained = true;
            
            f(m);
            A_CHECK("lGaps == 0 && rGaps == 0", "No gaps expected for a TP");
            
            stats.data[align.cID].tp++;
            x.aLvl.r2r[m->id()]++;
        }
        else
        {
            // At the alignment level, anything but a perfect match is a FP
            isTP = false;
            
            x.afp.push_back(align.name);
            
            // Can we at least match by overlapping?
            const auto m = stats.inters[align.cID].overlap(l);
            
            if (m)
            {
                f(m);
                A_ASSERT(lGaps != 0 || rGaps != 0);
                
                // Gap to the left?
                if (l.start < m->l().start)
                {
                    const auto gap = Locus(l.start, m->l().start-1);
                    
                    x.bLvl.fp->map(gap);
                }
                
                // Gap to the right?
                if (l.end > m->l().end)
                {
                    const auto gap = Locus(m->l().end+1, l.end);
                    
                    x.bLvl.fp->map(gap);
                }
                
                stats.data[align.cID].fp++;
            }
            else if (isMetaQuin(align.cID))
            {
                stats.data[align.cID].fp++;
            }
        }
    }
    
    if (isTP)
    {
        x.aLvl.m.tp()++;
    }
    else
    {
        x.aLvl.m.fp()++;
    }
}

MAlign::Stats MAlign::analyze(const FileName &file, const Options &o)
{
    auto stats = init();
    
    auto classify = [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        // Intron? Probably a mistake.
        if (info.skip)
        {
            o.warn("Skipped alignment: " + x.name);
        }
        
        if (!x.mapped)
        {
            return;
        }
        
        stats.update(x, isMetaQuin);
        
        if (info.skip)
        {
            return;
        }
        
        if (isMetaQuin(x.cID))
        {
            classifyAlign(stats, x);
        }
        else if (Standard::isGenomic(x.cID))
        {
            classifyAlign(stats, x);
        }
        else
        {
            o.warn(x.cID);
        }
    };
    
    o.analyze(file);
    
    ParserBAM::parse(file, [&](ParserBAM::Data &align, const ParserBAM::Info &info)
    {
        classify(align, info);
    });
    
    /*
     * -------------------- Calculating statistics --------------------
     */
    
    Base stp = 0;
    Base sfp = 0;
    Base gtp = 0;
    Base gfp = 0;
    
    // For each genome...
    for (const auto &i : stats.inters)
    {
        const auto &cID = i.first;
        
        // For each region...
        for (const auto &j : i.second.data())
        {
            const auto &rID  = j.first;
            const auto isSyn = isMetaQuin(cID);
            
            const auto m = stats.inters.at(cID).find(rID);
            A_ASSERT(m);
            
            // Statistics for the region
            const auto rs = m->stats();
            
            A_CHECK(rs.length, "Empty region: " + rID);
            
            if (isSyn)
            {
                stats.s2l[rID] = rs.length;
                stats.s2c[rID] = rs.nonZeros;
                
                // Sensitivty for the region
                stats.g2s[rID] = static_cast<Proportion>(stats.s2c.at(rID)) / stats.s2l.at(rID);
            }
            else
            {
                stats.g2l[rID] = rs.length;
                stats.g2c[rID] = rs.nonZeros;
                
                // Sensitivty for the gene
                stats.g2s[rID] = static_cast<Proportion>(stats.g2c.at(rID)) / stats.g2l.at(rID);
            }
            
            A_ASSERT(stats.s2c[rID] <= stats.s2l[rID]);
            A_ASSERT(stats.g2c[rID] <= stats.g2l[rID]);
            
            if (!stats.data.count(cID))
            {
                o.warn("No alignments found for " + cID);
            }
            else
            {
                // TP at the base level
                const auto btp = stats.data.at(cID).align.count(rID) ? stats.data.at(cID).align.at(rID) : 0;
                
                // FP at the base level (requires overlapping)
                const auto bfp = (stats.data.at(i.first).lGaps.count(rID) ? stats.data.at(i.first).lGaps.at(rID) : 0)
                +
                (stats.data.at(i.first).rGaps.count(rID) ? stats.data.at(i.first).rGaps.at(rID) : 0);
                
                A_ASSERT(!isnan(btp) && btp >= 0);
                A_ASSERT(!isnan(bfp) && bfp >= 0);
                
                // Precision at the base level
                const auto bpc = static_cast<Proportion>(btp) / (btp + bfp);
                
                A_ASSERT(isnan(bpc) || (bpc >= 0.0 && bpc <= 1.0));
                
                if (isMetaQuin(cID))
                {
                    stp += btp;
                    sfp += bfp;
                }
                else
                {
                    gtp += btp;
                    gfp += bfp;
                }
                
                stats.g2p[rID] = bpc;
            }
            
            A_ASSERT(stp >= 0);
            A_ASSERT(sfp >= 0);
        }
        
        auto &x = stats.data.at(cID);
        
        /*
         * Aggregating alignment statistics for the whole chromosome
         */
        
        const auto atp = x.aLvl.m.tp();
        const auto afp = x.aLvl.m.fp();
        
        if (isMetaQuin(cID))
        {
            stats.sa.tp() += atp;
            stats.sa.fp() += afp;
        }
        else
        {
            stats.ga.tp() += atp;
            stats.ga.fp() += afp;
        }
        
        /*
         * Aggregating base statistics for the whole chromosome
         */
        
        // Statistics for the reference region (TP)
        const auto ts = stats.inters.at(cID).stats();
        
        // Statistics for the non-reference region (FP)
        const auto fs = x.bLvl.fp->stats();
        
        if (isMetaQuin(cID))
        {
            stats.sb.tp() += ts.nonZeros;
            stats.sb.fp() += fs.nonZeros;
            stats.sb.fn() += ts.length - ts.nonZeros;
        }
        else
        {
            stats.gb.tp() += ts.nonZeros;
            stats.gb.fp() += fs.nonZeros;
            stats.gb.fn() += ts.length - ts.nonZeros;
        }
    }
    
    A_ASSERT(!stats.g2s.empty());
    A_ASSERT(!stats.s2l.empty());
    A_ASSERT(!stats.s2c.empty());
    
    A_ASSERT(stats.s2l.size() == stats.s2c.size());
    
    A_ASSERT(stats.sb.pc() >= 0.0 && stats.sb.pc() <= 1.0);
    A_ASSERT(isnan(stats.gb.pc()) || (stats.gb.pc() >= 0.0 && stats.gb.pc() <= 1.0));
    
    A_ASSERT(stats.sb.sn() >= 0.0 && stats.sb.sn() <= 1.0);
    A_ASSERT(isnan(stats.gb.sn()) || (stats.gb.sn() >= 0.0 && stats.gb.sn() <= 1.0));
    
    return stats;
}

static void writeSummary(const FileName &file,
                         const FileName &src,
                         const MAlign::Stats &stats,
                         const MAlign::Options &o)
{
    extern FileName BedRef();
    
    const auto &r = Standard::instance().r_meta;
    
    const auto sums2c = sum(stats.s2c);
    const auto sums2l = sum(stats.s2l);
    const auto sumg2c = sum(stats.g2c);
    const auto sumg2l = sum(stats.g2l);
    
    assert(sums2l >= sums2c);
    assert(sumg2l >= sumg2c);
    
    const auto summary = "-------MetaAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Alignment file: %2%\n\n"
                         "-------Alignments\n\n"
                         "       Synthetic: %3% (%4%)\n"
                         "       Genome:    %5% (%6%)\n"
                         "       Dilution:  %7$.2f\n\n"
                         "-------Reference annotation (Synthetic)\n\n"
                         "       Synthetic: %8% regions\n"
                         "       Synthetic: %9% bases\n\n"
                         "-------Reference annotation (Genome)\n\n"
                         "       Genome: %10% regions\n"
                         "       Genome: %11% bases\n\n"
                         "-------Comparison of alignments to annotation (Synthetic)\n\n"
                         "       *Alignment level\n"
                         "       Correct:     %12%\n"
                         "       Incorrect:   %13%\n\n"
                         "       Precision:   %14$.4f\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %15%\n"
                         "       Uncovered:   %16%\n"
                         "       Erroneous:   %17%\n"
                         "       Total:       %18%\n\n"
                         "       Sensitivity: %19$.4f\n"
                         "       Precision:   %20$.4f\n\n"
                         "-------Comparison of alignments to annotation (Genome)\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %21%\n"
                         "       Uncovered:   %22%\n"
                         "       Total:       %23%\n\n"
                         "       Sensitivity: %24$.4f\n";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()                // 1
                                            % src                     // 2
                                            % stats.nSeqs             // 3
                                            % (100 * stats.pSyn()) // 4
                                            % stats.nEndo             // 5
                                            % (100 * stats.pEndo()) // 6
                                            % stats.dilution()        // 7
                                            % r.nMicroSyn()           // 8
                                            % r.nBaseSyn()            // 9
                                            % r.nMicroGen()           // 10
                                            % r.nBaseGen()            // 11
                                            % stats.sa.tp()           // 12
                                            % stats.sa.fp()           // 13
                                            % stats.sa.pc()           // 14
                                            % stats.sb.tp()           // 15
                                            % stats.sb.fn()           // 16
                                            % stats.sb.fp()           // 17
                                            % (stats.sb.tp() + stats.sb.fp() + stats.sb.fn()) // 18
                                            % toString(stats.sb.sn())         // 19
                                            % toString(stats.sb.pc())         // 20
                                            % stats.gb.tp()                   // 21
                                            % stats.gb.fn()                   // 22
                                            % (stats.gb.tp() + stats.gb.fn()) // 23
                                            % toString(stats.gb.sn())         // 24
                     ).str());
    o.writer->close();
}

static void writeBQuins(const FileName &file,
                        const MAlign::Stats &stats,
                        const MAlign::Options &o)
{
#ifdef DEBUG_MALIGN
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ChrID" % "Position" % "Label").str());
    
    // For each chromosome...
    for (const auto &i : stats.data)
    {
        const auto &cID = i.first;
        
        // For each region in the chromosome...
        for (const auto &j : stats.inters.at(cID).data())
        {
            // For each mapped fragment in the region...
            for (const auto &k : j.second._data)
            {
                const auto pos = (toString(k.second.start) + "-" + toString(k.second.end));
                o.writer->write((boost::format(format) % cID % pos % "TP").str());
            }
            
            const auto zeros = j.second.zeros();
            
            for (const auto &k : zeros)
            {
                const auto pos = (toString(k.start) + "-" + toString(k.end));
                o.writer->write((boost::format(format) % cID % pos % "FN").str());
            }
        }
    }
    
    o.writer->close();
#endif
}

static void writeQuins(const FileName &file, const MAlign::Stats &stats, const MAlign::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4$.4f\t%5$.4f";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Reads"
                                           % "Sn"
                                           % "Pc").str());
    
    // For each chromosome...
    for (const auto &i : stats.inters)
    {
        const auto &cID = i.first;
        
        // Only the synthetic chromosome...
        if (isMetaQuin(cID))
        {
            // For each sequin region...
            for (const auto &j : i.second.data())
            {
                const auto &sID = j.first;
                
                // Data for the chromosome
                const auto &x = stats.data.at(cID);
                
                // Number of reads mapped to the region
                const auto reads = x.aLvl.r2r.count(sID) ? x.aLvl.r2r.at(sID) : 0;
                
                o.writer->write((boost::format(format) % sID
                                                       % stats.s2l.at(sID)
                                                       % reads
                                                       % stats.g2s.at(sID)
                                                       % stats.g2p.at(sID)).str());
            }
        }
    }
    
    o.writer->close();
}

static void writeQueries(const FileName &file, const MAlign::Stats &stats, const MAlign::Options &o)
{
#ifdef DEBUG_MALIGN
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "Reads" % "Label").str());
    
    for (const auto &i : stats.data)
    {
        if (isMetaQuin(i.first))
        {
            for (const auto &j : i.second.afp)
            {
                o.writer->write((boost::format(format) % j % "FP").str());
            }
        }
    }
    
    o.writer->close();
#endif
}

void MAlign::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = analyze(files, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating MetaAlign_summary.stats
     */
    
    writeSummary("MetaAlign_summary.stats", files[0], stats[0], o);
    
    /*
     * Generating MetaAlign_quins.stats
     */

    writeQuins("MetaAlign_sequins.csv", stats[0], o);
    
    /*
     * Generating MetaAlign_queries.stats (for debugging)
     */

    writeQueries("MetaAlign_queries.stats", stats[0], o);
    
    /*
     * Generating MetaAlign_rbase.stats (for debugging)
     */

    writeBQuins("MetaAlign_rbase.stats", stats[0], o);
}
