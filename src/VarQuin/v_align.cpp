#include "VarQuin/VarQuin.hpp"
#include "VarQuin/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

#ifdef DEBUG_VALIGN
static std::ofstream __bWriter__;
#endif

static void writeBase(const ChrID &cID, const Locus &l, const Label &label)
{
#ifdef DEBUG_VALIGN
    __bWriter__ << cID << "\t" << l.start << "\t" << l.end << "\t" << label << "\n";
#endif
}

static VAlign::Stats init()
{
    const auto &r = Standard::instance().r_var;

    VAlign::Stats stats;
    
    stats.inters = r.mInters();
    A_ASSERT(!stats.inters.empty());

    for (const auto &i : stats.inters)
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
        
        // Does the read aligned within a region?
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
            else if (isVarQuin(align.cID))
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
        x.aLvl.m.tp()++;
    }
    else
    {
        x.aLvl.m.fp()++;
    }
}

VAlign::Stats VAlign::analyze(const FileName &gen, const FileName &seqs, const Options &o)
{
    auto stats = init();
    
    o.info(std::to_string(stats.inters.size()) + " chromosomes in the reference");

#ifdef DEBUG_VALIGN
    __bWriter__.open(o.work + "/VarAlign_qbase.stats");
#endif

    auto classify = [&](ParserSAM::Data &x, const ParserSAM::Info &info)
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
        
        stats.update(x, isVarQuin);
        
        if (info.skip)
        {
            return;
        }
        
        if (isReverseGenome(x.cID))
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

    /*
     * Analyzing genomic alignments
     */
    
    o.analyze(gen);
    
    ParserSAM::parse(gen, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        if (!isReverseGenome(x.cID))
        {
            classify(x, info);
        }
    });
    
    /*
     * Analyzing sequin alignments (also in the forward genome)
     */
    
    o.analyze(seqs);
    
    ParserSAM::parse(seqs, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        /*
         * Important: the aligments will be on the forward genome. We must convert them to
         *            the reverse genome.
         */

        // Eg: chr2 to chrev2
        x.cID = toReverse(x.cID);
        
        if (isReverseGenome(x.cID))
        {
            classify(x, info);
        }
        else
        {
            o.logInfo("Invalid chromosome for sequins: " + x.cID + "." + x.name);
        }
    });

#ifdef DEBUG_VALIGN
    __bWriter__.close();
#endif

    o.info("Alignments analyzed. Generating statistics...");
    
    /*
     * -------------------- Calculating statistics --------------------
     */

    Base stp = 0;
    Base sfp = 0;
    Base gtp = 0;
    Base gfp = 0;
    
    o.info("Analyzing " + std::to_string(stats.inters.size()) + " chromsomes");

    // For each chromosome...
    for (const auto &i : stats.inters)
    {
        const auto &cID = i.first;
        
        // For each region...
        for (const auto &j : i.second.data())
        {
            const auto &rID  = j.first;
            const auto isSyn = isReverseGenome(cID);

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
                
                if (isSyn)
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
        
        if (isReverseGenome(cID))
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
        
        if (isReverseGenome(cID))
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

    //A_ASSERT(!stats.g2r.empty());
    A_ASSERT(!stats.g2s.empty());
    A_ASSERT(!stats.s2l.empty());
    A_ASSERT(!stats.s2c.empty());

    A_ASSERT(stats.s2l.size() == stats.s2c.size());
    //A_ASSERT(stats.g2r.size() == stats.g2s.size());

    A_ASSERT(stats.sb.pc() >= 0.0 && stats.sb.pc() <= 1.0);
    A_ASSERT(isnan(stats.gb.pc()) || (stats.gb.pc() >= 0.0 && stats.gb.pc() <= 1.0));

    A_ASSERT(stats.sb.sn() >= 0.0 && stats.sb.sn() <= 1.0);
    A_ASSERT(isnan(stats.gb.sn()) || (stats.gb.sn() >= 0.0 && stats.gb.sn() <= 1.0));
    
    return stats;
}

void VAlign::writeSummary(const FileName &file,
                          const FileName &gen,
                          const FileName &seq,
                          const VAlign::Stats &stats,
                          const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;

    const auto sums2c = sum(stats.s2c);
    const auto sums2l = sum(stats.s2l);
    const auto sumg2c = sum(stats.g2c);
    const auto sumg2l = sum(stats.g2l);
    
    A_ASSERT(sums2l >= sums2c);
    A_ASSERT(sumg2l >= sumg2c);

    const auto summary = "-------VarAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Genome alignment file: %2%\n"
                         "       Synthetic alignment file: %3%\n\n"
                         "-------Alignments\n\n"
                         "       Synthetic: %4% (%5%%%)\n"
                         "       Genome:    %6% (%7%%%)\n"
                         "       Dilution:  %8$.4f%%\n\n"
                         "-------Reference annotation (Synthetic)\n\n"
                         "       Synthetic: %9% regions\n"
                         "       Synthetic: %10% bases\n\n"
                         "-------Reference annotation (Genome)\n\n"
                         "       Genome: %11% regions\n"
                         "       Genome: %12% bases\n\n"
                         "-------Comparison of alignments to annotation (Synthetic)\n\n"
                         "       *Alignment level\n"
                         "       Inside regions:  %13%\n"
                         "       Outside regions: %14%\n\n"
                         "       Precision:      %15$.4f\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %16%\n"
                         "       Uncovered:   %17%\n"
                         "       Erroneous:   %18%\n"
                         "       Total:       %19%\n\n"
                         "       Sensitivity: %20$.4f\n"
                         "       Precision:   %21$.4f\n\n"
                         "-------Comparison of alignments to annotation (Genome)\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %22%\n"
                         "       Uncovered:   %23%\n"
                         "       Total:       %24%\n\n"
                         "       Sensitivity: %25$.4f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % o.rBed           // 1
                                            % gen              // 2
                                            % seq              // 3
                                            % stats.nSyn       // 4
                                            % stats.propSyn()  // 5
                                            % stats.nGen       // 6
                                            % stats.propGen()  // 7
                                            % (100.0 * stats.dilution()) // 8
                                            % r.nGeneSyn()     // 9
                                            % r.nBaseSyn()     // 10
                                            % r.nGeneGen()     // 11
                                            % r.nBaseGen()     // 12
                                            % stats.sa.tp()    // 13 (inside region)
                                            % stats.sa.fp()    // 14 (outside region)
                                            % stats.sa.pc()    // 15
                                            % stats.sb.tp()    // 16
                                            % stats.sb.fn()    // 17
                                            % stats.sb.fp()    // 18
                                            % (stats.sb.tp() + stats.sb.fp() + stats.sb.fn()) // 19
                                            % stats.sb.sn()    // 20
                                            % stats.sb.pc()    // 21
                                            % stats.gb.tp()    // 22
                                            % stats.gb.fn()    // 23
                                            % (stats.gb.tp() + stats.gb.fn()) // 24
                                            % stats.gb.sn()    // 25
                     ).str());
    o.writer->close();
    
    // Only for generating a report, easier for parsing. (undocumented)
    if (o.report)
    {
        #define WRITE(x,y) o.writer->write((boost::format("%1%\t%2%")   % x % y).str());

        o.generate("VarAlign_report.csv");
        o.writer->open("VarAlign_report.csv");
        
        WRITE("ASyn", stats.nSyn);
        WRITE("ASynP", stats.propSyn());
        WRITE("AGen", stats.nGen);
        WRITE("AGenP", stats.propGen());
        WRITE("Dilution", (100.0 * stats.dilution()));

        WRITE("nGeneSyn", r.nGeneSyn());
        WRITE("nBaseSyn", r.nBaseSyn());
        
        WRITE("RegionTP", stats.sa.tp());
        WRITE("RegionFP", stats.sa.fp());
        WRITE("RegionPC", stats.sa.pc());
        
        WRITE("BaseTP", stats.sb.tp());
        WRITE("BaseFN", stats.sb.fn());
        WRITE("BaseFP", stats.sb.fp());
        WRITE("BaseSN", stats.sb.sn());
        WRITE("BasePC", stats.sb.pc());
        
        WRITE("GenomeTP", stats.gb.tp());
        WRITE("GenomeFN", stats.gb.fn());
        WRITE("GenomeSN", stats.gb.sn());
        
        o.writer->close();
    }
}

void VAlign::writeBQuins(const FileName &file,
                         const VAlign::Stats &stats,
                         const VAlign::Options &o)
{
#ifdef DEBUG_VALIGN
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

void VAlign::writeQuins(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4$.4f\t%5$.4f";
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Reads"
                                           % "Sn"
                                           % "Pc").str());

    o.logInfo("writeQuins: " + std::to_string(stats.inters.size()));
    
    // For each chromosome...
    for (const auto &i : stats.inters)
    {
        const auto &cID = i.first;
        
        // Only the synthetic chromosome...
        if (isVarQuin(cID))
        {
            o.logInfo(i.first + " - " + std::to_string(i.second.data().size()));
            
            // For each sequin region...
            for (const auto &j : i.second.data())
            {
                o.logInfo(j.first);
                
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

void VAlign::writeQueries(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
#ifdef DEBUG_VALIGN
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "Reads" % "Label").str());

    for (const auto &i : stats.data)
    {
        if (isVarQuin(i.first))
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

void VAlign::report(const FileName &gen, const FileName &seqs, const Options &o)
{
    const auto stats = analyze(gen, seqs, o);

    o.info("Generating statistics");
    
    /*
     * Generating VarAlign_summary.stats
     */
    
    VAlign::writeSummary("VarAlign_summary.stats", gen, seqs, stats, o);

    /*
     * Generating VarAlign_sequins.csv
     */
    
    if (!o.report)
    {
        VAlign::writeQuins("VarAlign_sequins.csv", stats, o);
    }

    /*
     * Generating VarAlign_queries.stats (for debugging)
     */
    
    VAlign::writeQueries("VarAlign_queries.stats", stats, o);

    /*
     * Generating VarAlign_rbase.stats (for debugging)
     */
    
    VAlign::writeBQuins("VarAlign_rbase.stats", stats, o);
}
