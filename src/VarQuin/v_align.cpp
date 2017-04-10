#include "tools/tools.hpp"
#include "VarQuin/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

extern FileName BedRef();

#ifdef DEBUG_VALIGN
static std::ofstream __bWriter__;
#endif

static void writeBase(const ChrID &cID, const Locus &l, const Label &label)
{
#ifdef DEBUG_VALIGN
    __bWriter__ << cID << "\t" << l.start << "\t" << l.end << "\t" << label << "\n";
#endif
}

static void classifyAlign(VAlign::Performance &stats, ParserSAM::Data &align)
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

            stats.data[align.cID].fp++;
                
            // Not overlapping with the reference regions
            writeBase(align.cID, l, "FP");
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

VAlign::Stats VAlign::analyze(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto &r = Standard::instance().r_var;

    VAlign::Stats stats;
    
    auto initPerf = [&](Performance &p)
    {
        p.inters = r.mInters();
        A_ASSERT(!p.inters.empty());
        
        for (const auto &i : p.inters)
        {
            const auto &cID = i.first;
            
            /*
             * We'd like to know the length of the chromosome but we don't have the information.
             * It doesn't matter because we can simply initalize it to the the maximum possible.
             * The data structure in MergedInterval will be efficient not to waste memory.
             */
            
            MergedInterval *mi = new MergedInterval(cID, Locus(1, std::numeric_limits<Base>::max()));
            p.data[cID].bLvl.fp = std::shared_ptr<MergedInterval>(mi);
        }
        
        return stats;
    };
    
    stats.endo = std::shared_ptr<Performance>(new Performance());
    stats.seqs = std::shared_ptr<Performance>(new Performance());
    
    initPerf(*(stats.endo));
    initPerf(*(stats.seqs));
    
    o.logInfo(std::to_string(stats.seqs->inters.size()) + " regions");

#ifdef DEBUG_VALIGN
    __bWriter__.open(o.work + "/VarAlign_qbase.stats");
#endif

    auto classify = [&](ParserSAM::Data &x, const ParserSAM::Info &info, Performance &p)
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
            p.nNA++;
            return;
        }
        
        p.nMap++;
        
        if (info.skip)
        {
            return;
        }
        
        classifyAlign(p, x);
    };

    /*
     * Analyzing endogenous alignments
     */
    
    o.analyze(endo);
    
    ParserSAM::parse(endo, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        classify(x, info, *(stats.endo));
    });
    
    /*
     * Analyzing sequin alignments
     */
    
    o.analyze(seqs);
    
    ParserSAM::parse(seqs, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
    {
        classify(x, info, *(stats.seqs));
    });

#ifdef DEBUG_VALIGN
    __bWriter__.close();
#endif

    o.info("Alignments analyzed. Generating statistics...");
    
    /*
     * -------------------- Calculating statistics --------------------
     */

    auto analyze = [&](Performance &p)
    {
        Base tp = 0;
        Base fp = 0;
        
        o.info("Analyzing " + std::to_string(p.inters.size()) + " regions");
        
        // For each chromosome...
        for (const auto &i : p.inters)
        {
            const auto &cID = i.first;
            
            // For each reference region...
            for (const auto &j : i.second.data())
            {
                const auto &rID = j.first;
                
                const auto m = p.inters.at(cID).find(rID);
                A_ASSERT(m);
                
                // Statistics for the region
                const auto rs = m->stats();
                
                A_CHECK(rs.length, "Empty region: " + rID);
                
                p.length[rID] = rs.length;
                p.covered[rID] = rs.nonZeros;
                
                // Sensitivty for the region
                p.r2s[rID] = static_cast<Proportion>(p.covered.at(rID)) / p.length.at(rID);
                
                A_ASSERT(p.covered[rID] <= p.length[rID]);
                
                if (!p.data.count(cID))
                {
                    o.warn("No alignments found for: " + cID);
                }
                else
                {
                    // TP at the base level
                    const auto btp = p.data.at(cID).align.count(rID) ? p.data.at(cID).align.at(rID) : 0;
                    
                    // FP at the base level (requires overlapping)
                    const auto bfp = (p.data.at(i.first).lGaps.count(rID) ? p.data.at(i.first).lGaps.at(rID) : 0)
                                                                          +
                                     (p.data.at(i.first).rGaps.count(rID) ? p.data.at(i.first).rGaps.at(rID) : 0);
                    
                    A_ASSERT(!isnan(btp) && btp >= 0);
                    A_ASSERT(!isnan(bfp) && bfp >= 0);
                    
                    // Precision at the base level
                    const auto bpc = static_cast<Proportion>(btp) / (btp + bfp);
                    
                    A_ASSERT(isnan(bpc) || (bpc >= 0.0 && bpc <= 1.0));
                    
                    tp += btp;
                    fp += bfp;
                    
                    p.r2p[rID] = bpc;
                }
                
                A_ASSERT(tp >= 0);
                A_ASSERT(fp >= 0);
            }
            
            auto &x = p.data.at(cID);
            
            /*
             * Aggregating alignment statistics for the whole chromosome
             */
            
            const auto atp = x.aLvl.m.tp();
            const auto afp = x.aLvl.m.fp();
            
            p.align.tp() += atp;
            p.align.fp() += afp;
            
            /*
             * Aggregating base statistics for the whole chromosome
             */
            
            // Statistics for the reference region (TP)
            const auto ts = p.inters.at(cID).stats();
            
            // Statistics for the non-reference region (FP)
            const auto fs = x.bLvl.fp->stats();
            
            p.base.tp() += ts.nonZeros;
            p.base.fp() += fs.nonZeros;
            p.base.fn() += ts.length - ts.nonZeros;
        }
        
        //A_ASSERT(!p.g2r.empty());
        A_ASSERT(!p.r2s.empty());
        A_ASSERT(!p.length.empty());
        A_ASSERT(!p.covered.empty());
        
        A_ASSERT(p.length.size() == p.covered.size());
        //A_ASSERT(p.g2r.size() == p.g2s.size());
    };

    analyze(*(stats.seqs));
    analyze(*(stats.endo));
    
    const auto x = stats.endo;
    const auto y = stats.seqs;

    A_ASSERT(isnan(x->base.pc()) || (x->base.pc() >= 0.0 && x->base.pc() <= 1.0));
    A_ASSERT(isnan(x->base.sn()) || (x->base.sn() >= 0.0 && x->base.sn() <= 1.0));

    A_ASSERT(y->base.pc() >= 0.0 && y->base.pc() <= 1.0);
    A_ASSERT(y->base.sn() >= 0.0 && y->base.sn() <= 1.0);

    return stats;
}

void VAlign::writeSummary(const FileName &file,
                          const FileName &gen,
                          const FileName &seq,
                          const VAlign::Stats &stats,
                          const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;

    const auto sumg2c = sum(stats.endo->covered);
    const auto sumg2l = sum(stats.endo->length);
    const auto sums2c = sum(stats.seqs->covered);
    const auto sums2l = sum(stats.seqs->length);
    
    A_ASSERT(sums2l >= sums2c);
    A_ASSERT(sumg2l >= sumg2c);
    
    const auto &endo = stats.endo;
    const auto &seqs = stats.seqs;

    const auto summary = "-------VarAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       Genome alignment file: %2%\n"
                         "       Synthetic alignment file: %3%\n\n"
                         "-------Alignments\n\n"
                         "       Genome:    %4% (%5%%%)\n"
                         "       Synthetic: %6% (%7%%%)\n"
                         "       Dilution:  %8$.4f%%\n\n"
                         "-------Reference regions\n\n"
                         "       Regions: %9% regions\n"
                         "       Regions: %10% bases\n\n"
                         "-------Comparison of alignments to annotation (Synthetic)\n\n"
                         "       *Alignment level\n"
                         "       Inside regions:  %11%\n"
                         "       Outside regions: %12%\n\n"
                         "       Precision:      %13$.4f\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %14%\n"
                         "       Uncovered:   %15%\n"
                         "       Erroneous:   %16%\n"
                         "       Total:       %17%\n\n"
                         "       Sensitivity: %18$.4f\n"
                         "       Precision:   %19$.4f\n\n"
                         "-------Comparison of alignments to annotation (Genome)\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %20%\n"
                         "       Uncovered:   %21%\n"
                         "       Total:       %22%\n\n"
                         "       Sensitivity: %23$.4f\n";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()              // 1
                                            % gen                   // 2
                                            % seq                   // 3
                                            % endo->nMap            // 4
                                            % stats.pEndo()         // 5
                                            % seqs->nMap            // 6
                                            % stats.pSeqs()         // 7
                                            % (100 * stats.pSeqs()) // 8
                                            % r.nRegs()             // 9
                                            % r.lRegs()             // 10
                                            % seqs->align.tp()      // 11
                                            % seqs->align.fp()      // 12
                                            % seqs->align.pc()      // 13
                                            % seqs->base.tp()       // 14
                                            % seqs->base.fn()       // 15
                                            % seqs->base.fp()       // 16
                                            % (seqs->base.tp() + seqs->base.fp() + seqs->base.fn()) // 17
                                            % seqs->base.sn()       // 18
                                            % seqs->base.pc()       // 19
                                            % endo->base.tp()       // 20
                                            % endo->base.fn()       // 21
                                            % (endo->base.tp() + endo->base.fn()) // 22
                                            % endo->base.sn()       // 23
                     ).str());
    o.writer->close();
    
    // Only for generating a report, easier for parsing. (undocumented)
    if (o.report)
    {
        #define WRITE(x,y) o.writer->write((boost::format("%1%\t%2%")   % x % y).str());

        o.generate("VarAlign_report.csv");
        o.writer->open("VarAlign_report.csv");
        
//        WRITE("ASyn", seqs.nSeqs);
//        WRITE("ASynP", seqs.propSyn());
//        WRITE("AGen", seqs.nEndo);
//        WRITE("AGenP", seqs.propGen());
//        WRITE("Dilution", (100.0 * seqs.dilution()));

        WRITE("nRegs", r.nRegs());
        WRITE("lRegs", r.lRegs());
        
        WRITE("RegionTP", seqs->align.tp());
        WRITE("RegionFP", seqs->align.fp());
        WRITE("RegionPC", seqs->align.pc());
        
        WRITE("BaseTP", seqs->base.tp());
        WRITE("BaseFN", seqs->base.fn());
        WRITE("BaseFP", seqs->base.fp());
        WRITE("BaseSN", seqs->base.sn());
        WRITE("BasePC", seqs->base.pc());
        
        WRITE("GenomeTP", endo->base.tp());
        WRITE("GenomeFN", endo->base.fn());
        WRITE("GenomeSN", endo->base.sn());
        
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

    o.logInfo("writeQuins: " + std::to_string(stats.seqs->inters.size()));
    
    // For each chromosome...
    for (const auto &i : stats.seqs->inters)
    {
        const auto &cID = i.first;
        
        o.logInfo(i.first + " - " + std::to_string(i.second.data().size()));
        
        // For each sequin region...
        for (const auto &j : i.second.data())
        {
            o.logInfo(j.first);
            
            const auto &sID = j.first;
            
            // Data for the chromosome
            const auto &x = stats.seqs->data.at(cID);
            
            // Number of reads mapped to the region
            const auto reads = x.aLvl.r2r.count(sID) ? x.aLvl.r2r.at(sID) : 0;
            
            o.writer->write((boost::format(format) % sID
                                                   % stats.seqs->length.at(sID)
                                                   % reads
                                                   % stats.seqs->r2s.at(sID)
                                                   % stats.seqs->r2p.at(sID)).str());
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
        for (const auto &j : i.second.afp)
        {
            o.writer->write((boost::format(format) % j % "FP").str());
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
    
    o.logInfo("VarAlign completed");
}
