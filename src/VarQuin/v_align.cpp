#include "VarQuin/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

static std::ofstream __bWriter__;

static void writeBase(const ChrID &cID, const Locus &l, const Label &label)
{
    __bWriter__ << cID << "\t" << l.start << "\t" << l.end << "\t" << label << "\n";
}

static VAlign::Stats init()
{
    const auto &r = Standard::instance().r_var;

    VAlign::Stats stats;
    
    stats.hist   = r.hist();
    stats.inters = r.mInters();

    std::cout << "[INFO]: " << stats.inters.size() << " chromosomes in the reference" << std::endl;
    
    assert(!stats.hist.empty());
    assert(!stats.inters.empty());

    for (const auto &i : stats.hist)
    {
        const auto &cID = i.first;
        
        /*
         * We'd like to know the length of the chromosome but we don't have the information.
         * It doesn't matter because we can simply initalize it to the the maximum possible.
         * The data structure in MergedInterval will be efficient not to waste memory.
         */
        
        MergedInterval *mi = new MergedInterval(cID, Locus(1, std::numeric_limits<Base>::max()));
        stats.data[cID].bLvl.fp = std::shared_ptr<MergedInterval>(mi);
    }

    return stats;
}

static void classifyAlign(VAlign::Stats &stats, ParserSAM::Data &align)
{
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
            
            assert(covered >= 0);
            assert(l.length() > lGaps);
            assert(l.length() > rGaps);
        };
        
        // Does the read aligned within a gene (or a region)?
        const auto m = stats.inters.at(align.cID).contains(l);
        
        if (m)
        {
            isContained = true;
            
            f(m);
            assert(lGaps == 0 && rGaps == 0);
            
            stats.data[align.cID].tp++;
            //stats.hist.at(align.cID).at(m->name())++;
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
                assert(lGaps != 0 || rGaps != 0);
                
                // Gap to the left?
                if (l.start < m->l().start)
                {
                    const auto gap = Locus(l.start, m->l().start-1);
                    
                    x.bLvl.fp->map(gap);
                    writeBase(align.cID, gap, "FP");
                }
                
                // Gap to the right?
                if (l.end > m->l().end)
                {
                    const auto gap = Locus(m->l().end+1, l.end);
                    
                    x.bLvl.fp->map(gap);
                    writeBase(align.cID, gap, "FP");
                }
                
                stats.data[align.cID].fp++;
                
                writeBase(align.cID, l, "FP");
            }
            else if (Standard::isSynthetic(align.cID))
            {
                stats.data[align.cID].fp++;
                
                /*
                 * The read is not aligned within the reference regions. We don't know whether this is
                 * a TP or FP.
                 */
                
                writeBase(align.cID, l, "FP");
            }
        }
    }
    
    if (isTP)
    {
        x.aLvl.tp()++;
    }
    else
    {
        x.aLvl.fp()++;
    }
}

VAlign::Stats VAlign::analyze(const FileName &file, const Options &o)
{
    auto stats = init();
    o.analyze(file);

    __bWriter__.open(o.work + "/VarAlign_qbase.stats");

    ParserSAM::parse(file, [&](ParserSAM::Data &align, const ParserSAM::Info &info)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        // Intron? Probably a mistake.
        if (info.skip)
        {
            o.warn("Skipped alignment: " + align.name);
        }
        
        if (!align.mapped)
        {
            return;
        }

        stats.update(align);

        if (info.skip)
        {
            return;
        }
        
        if (Standard::isSynthetic(align.cID))
        {
            //std::cout << k++ << std::endl;
            classifyAlign(stats, align);
        }
        else if (Standard::isGenomic(align.cID))
        {
            classifyAlign(stats, align);
        }
        else
        {
            o.warn(align.cID);
        }
    });
    
    __bWriter__.close();

    o.info("Alignments analyzed. Generating statistics.");
    
    /*
     * -------------------- Calculating statistics --------------------
     */

    Base stp = 0;
    Base sfp = 0;
    Base gtp = 0;
    Base gfp = 0;
    
    o.info("Analyzing " + toString(stats.hist.size()) + " chromsomes");

    // For each chromosome...
    for (const auto &i : stats.hist)
    {
        const auto &cID = i.first;
        
        auto &x = stats.data.at(cID);

        // For each region... (whole chromosome for the whole genome sequencing)
        for (const auto &j : i.second)
        {
            const auto &gID = j.first;

            // Reads aligned to the region
            stats.g2r[gID] = j.second;

            const auto isSyn = Standard::isSynthetic(cID);
            
            const auto m = stats.inters.at(cID).find(gID);
            assert(m);

            // Statistics for the gene (created by the interval)
            const auto ms = m->stats();

            assert(ms.length);
            
            if (isSyn)
            {
                stats.s2l[gID] = ms.length;
                stats.s2c[gID] = ms.nonZeros;

                // Sensitivty for the gene
                stats.g2s[gID] = static_cast<Proportion>(stats.s2c.at(gID)) / stats.s2l.at(gID);
            }
            else
            {
                stats.g2l[gID] = ms.length;
                stats.g2c[gID] = ms.nonZeros;

                // Sensitivty for the gene
                stats.g2s[gID] = static_cast<Proportion>(stats.g2c.at(gID)) / stats.g2l.at(gID);
            }
            
            assert(stats.s2c[gID] <= stats.s2l[gID]);
            assert(stats.g2c[gID] <= stats.g2l[gID]);

            if (!stats.data.count(cID))
            {
                o.warn("No alignments found for " + cID);
            }
            else
            {
                // TP at the base level
                const auto btp = stats.data.at(cID).align.count(gID) ? stats.data.at(cID).align.at(gID) : 0;
                
                // FP at the base level (requires overlapping)
                const auto bfp = (stats.data.at(i.first).lGaps.count(gID) ? stats.data.at(i.first).lGaps.at(gID) : 0)
                                                +
                                 (stats.data.at(i.first).rGaps.count(gID) ? stats.data.at(i.first).rGaps.at(gID) : 0);
                
                assert(!isnan(btp) && btp >= 0);
                assert(!isnan(bfp) && bfp >= 0);
                
                // Precison at the base level
                const auto bpc = static_cast<Proportion>(btp) / (btp + bfp);
                
                assert(isnan(bpc) || (bpc >= 0.0 && bpc <= 1.0));
                
                if (Standard::isSynthetic(cID))
                {
                    stp += btp;
                    sfp += bfp;
                }
                else
                {
                    gtp += btp;
                    gfp += bfp;
                }
                
                stats.g2p[gID] = bpc;
            }
            
            assert(stp >= 0);
            assert(sfp >= 0);
        }

        /*
         * Aggregating alignment statistics for the whole chromosome
         */
        
        const auto atp = x.aLvl.tp();
        const auto afp = x.aLvl.fp();
        
        if (Standard::isSynthetic(cID))
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
        
        if (Standard::isSynthetic(cID))
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

    assert(!stats.g2r.empty());
    assert(!stats.g2s.empty());
    assert(!stats.s2l.empty());
    assert(!stats.s2c.empty());

    assert(stats.s2l.size() == stats.s2c.size());
    assert(stats.g2r.size() == stats.g2s.size());

    assert(stats.sb.pc() >= 0.0 && stats.sb.pc() <= 1.0);
    assert(isnan(stats.gb.pc()) || (stats.gb.pc() >= 0.0 && stats.gb.pc() <= 1.0));

    assert(stats.sb.sn() >= 0.0 && stats.sb.sn() <= 1.0);
    assert(isnan(stats.gb.sn()) || (stats.gb.sn() >= 0.0 && stats.gb.sn() <= 1.0));
    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VAlign::Stats &stats, const VAlign::Options &o)
{
    extern FileName BedRef();

    const auto &r = Standard::instance().r_var;

    const auto sums2c = sum(stats.s2c);
    const auto sums2l = sum(stats.s2l);
    const auto sumg2c = sum(stats.g2c);
    const auto sumg2l = sum(stats.g2l);
    
    assert(sums2l >= sums2c);
    assert(sumg2l >= sumg2c);

    const auto summary = "-------VarAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       User alignment file: %2%\n\n"
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
                         "       Precison:    %14$.4f\n\n"
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
                                            % stats.n_syn             // 3
                                            % (100 * stats.synProp()) // 4
                                            % stats.n_gen             // 5
                                            % (100 * stats.genProp()) // 6
                                            % stats.dilution()        // 7
                                            % r.countGeneSyn()        // 8
                                            % r.countBaseSyn()        // 9
                                            % r.countGeneGen()        // 10
                                            % r.countBaseGen()        // 11
                                            % stats.sa.tp()           // 12
                                            % stats.sa.fp()           // 13
                                            % stats.sa.pc()           // 14
                                            % stats.sb.tp()                   // 15
                                            % stats.sb.fn()                   // 16
                                            % stats.sb.fp()                   // 17
                                            % (stats.sb.tp() + stats.sb.fp() + stats.sb.fn()) // 18
                                            % stats.sb.sn()                   // 19
                                            % stats.sb.pc()                   // 20
                                            % stats.gb.tp()                   // 21
                                            % stats.gb.fn()                   // 22
                                            % (stats.gb.tp() + stats.gb.fn()) // 23
                                            % stats.gb.sn()                   // 24
                     ).str());
    o.writer->close();
}

static void writeBQuins(const FileName &file,
                        const VAlign::Stats &stats,
                        const VAlign::Options &o)
{
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ChrID" % "Position" % "Label").str());
    
    for (const auto &i : stats.data)
    {
        const auto &cID = i.first;
        
        //if (Standard::isSynthetic(cID))
        {
            for (const auto &j : stats.inters.at(cID).data())
            {
                for (const auto &k : j.second._data)
                {
                    const auto pos = (toString(k.second.start) + "-" + toString(k.second.end));
                    
                    o.writer->write((boost::format(format) % cID
                                                           % pos
                                                           % "TP").str());
                }
                
                const auto zeros = j.second.zeros();
                
                for (const auto &k : zeros)
                {
                    const auto pos = (toString(k.start) + "-" + toString(k.end));
                    
                    o.writer->write((boost::format(format) % cID
                                                           % pos
                                                           % "FN").str());
                }
            }
        }
    }
    
    o.writer->close();
}

static void writeQuins(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Reads"
                                           % "Sn"
                                           % "Pc").str());

    // For each chromosome...
    for (const auto &i : stats.hist)
    {
        // Only the synthetic chromosome...
        if (Standard::isSynthetic(i.first))
        {
            for (const auto &j : i.second)
            {
                o.writer->write((boost::format(format) % j.first
                                                       % stats.s2l.at(j.first)
                                                       % stats.g2r.at(j.first)
                                                       % stats.g2s.at(j.first)
                                                       % stats.g2p.at(j.first)).str());
            }
        }
    }

    o.writer->close();
}

static void writeQueries(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "Reads" % "Label").str());

    for (const auto &i : stats.data)
    {
        if (Standard::isSynthetic(i.first))
        {
            for (const auto &j : i.second.afp)
            {
                o.writer->write((boost::format(format) % j % "FP").str());
            }
        }
    }
    
    o.writer->close();
}

void VAlign::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, o);

    o.info("Generating statistics");
    
    /*
     * Generating VarAlign_summary.stats
     */
    
    writeSummary("VarAlign_summary.stats", file, stats, o);

    /*
     * Generating VarAlign_quins.stats
     */
    
    writeQuins("VarAlign_sequins.csv", stats, o);

    /*
     * Generating VarAlign_queries.stats
     */
    
    writeQueries("VarAlign_queries.stats", stats, o);

    /*
     * Generating VarAlign_rbase.stats
     */
    
    writeBQuins("VarAlign_rbase.stats", stats, o);
    
    /*
     * Generating VarAlign_report.pdf
     */

    o.report->open("VarAlign_report.pdf");
    o.report->addTitle("VarAlign");
    o.report->addFile("VarAlign_summary.stats");
    o.report->addFile("VarAlign_sequins.csv");
    o.report->addFile("VarAlign_queries.stats");
}