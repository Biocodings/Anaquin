#include "tools/tools.hpp"
#include "VarQuin/v_align.hpp"
#include "parsers/parser_bam.hpp"
#include "writers/json_writer.hpp"

using namespace Anaquin;

extern FileName Bed1Ref();

#ifdef DEBUG_VALIGN
#include <fstream>
static std::ofstream __bWriter__;
#endif

static void writeBase(const ChrID &cID, const Locus &l, const Label &label)
{
#ifdef DEBUG_VALIGN
    __bWriter__ << cID << "\t" << l.start << "\t" << l.end << "\t" << label << "\n";
#endif
}

static void classifyAlign(VAlign::Performance &stats, ParserBAM::Data &align)
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
    
    auto initP = [&](Performance &p)
    {
        p.inters = r.mInters();
        A_ASSERT(!p.inters.empty());
        
        // For each chromosome...
        for (const auto &i : p.inters)
        {
            const auto &cID = i.first;
            
            /*
             * We'd like to know the length of the chromosome but we don't have the information.
             * It doesn't matter because we can initalize it to the the maximum possible.
             * The data structure in MergedInterval will be efficient not to waste memory.
             */
            
            p.data[cID].bLvl.fp = std::shared_ptr<MergedInterval>(
                    new MergedInterval(cID, Locus(1, std::numeric_limits<Base>::max())));
        }
        
        return stats;
    };
    
    stats.endo = std::shared_ptr<Performance>(new Performance());
    stats.seqs = std::shared_ptr<Performance>(new Performance());
    
    initP(*(stats.endo));
    initP(*(stats.seqs));
    
#ifdef DEBUG_VALIGN
    __bWriter__.open(o.work + "/VarAlign_qbase.stats");
#endif

    const auto r2 = r.regs2();
    
    auto classify = [&](ParserBAM::Data &x, const ParserBAM::Info &info, Performance &p)
    {
        if (info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        else if (!contains(r2, x.cID, x.l))
        {
            return;
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
    
    if (!endo.empty())
    {
        /*
         * Analyzing endogenous alignments
         */
        
        o.analyze(endo);
        
        ParserBAM::parse(endo, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
        {
            classify(x, info, *(stats.endo));
        });
    }

    /*
     * Analyzing sequin alignments
     */
    
    o.analyze(seqs);
    
    ParserBAM::parse(seqs, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
    {
        classify(x, info, *(stats.seqs));
    });

#ifdef DEBUG_VALIGN
    __bWriter__.close();
#endif

    /*
     * -------------------- Calculating statistics --------------------
     */

    auto analyze = [&](Performance &p)
    {
        Base tp = 0;
        Base fp = 0;
        
        o.info("Analyzing " + std::to_string(p.inters.size()) + " chromosomes");
        
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
                
                p.length[rID]  = rs.length;
                p.covered[rID] = rs.nonZeros;
                
                A_ASSERT(p.length.at(rID) >= 2.0 * o.edge);
                
                // Sensitivty for the region
                p.r2s[rID] = static_cast<Proportion>(p.covered.at(rID)) / (p.length.at(rID) - 2.0 * o.edge);
                
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
            p.base.fn() += (ts.length - (p.inters.at(cID).size() * 2.0 * o.edge)) - ts.nonZeros;
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

//    A_ASSERT(y->base.pc() >= 0.0 && y->base.pc() <= 1.0);
//    A_ASSERT(y->base.sn() >= 0.0 && y->base.sn() <= 1.0);

    return stats;
}

static std::map<std::string, std::string> jsonD(const FileName &endo,
                                                const FileName &seqs,
                                                const VAlign::Stats &stats,
                                                const VAlign::Options &o)
{
    const auto &r = Standard::instance().r_var;

    std::map<std::string, std::string> x;

    #define S0(x) (isnan(x) ? "-" : toString(x,0))
    #define S2(x) (isnan(x) ? "-" : toString(x,2))
    
    x["ref"]      = Bed1Ref();
    x["endo"]     = endo;
    x["seqs"]     = seqs;
    x["nNA"]      = stats.seqs->nNA;
    x["nEndo"]    = S0(stats.endo ? stats.endo->nMap : NAN);
    x["pEndo"]    = S2(stats.endo ? stats.pEndo() : NAN);
    x["nSeqs"]    = S0(stats.seqs->nMap);
    x["pSeqs"]    = S2(stats.pSeqs());
    x["dilution"] = S2(100 * stats.pSeqs());
    x["nRegs"]    = S0(r.nRegs());
    x["lRegs"]    = S0(r.lRegs());
    x["aTP"]      = S0(stats.seqs->align.tp());
    x["aFP"]      = S0(stats.seqs->align.fp());
    x["aPC"]      = S2(stats.seqs->align.pc());
    x["bTP"]      = S0(stats.seqs->base.tp());
    x["bFP"]      = S0(stats.seqs->base.fp());
    x["bFN"]      = S2(stats.seqs->base.fn());
    x["bPC"]      = S2(stats.seqs->base.pc());
    x["bSN"]      = S2(stats.seqs->base.sn());
    x["bPC"]      = S2(stats.seqs->base.pc());
    x["bTot"]     = S0(stats.seqs->base.tp() + stats.seqs->base.fp() + stats.seqs->base.fn());
    x["eTP"]      = S0(stats.endo ? stats.endo->base.tp() : NAN);
    x["eFN"]      = S0(stats.endo ? stats.endo->base.fn() : NAN);
    x["eSN"]      = S2(stats.endo ? stats.endo->base.sn() : NAN);
    x["eTot"]     = S0(stats.endo ? stats.endo->base.tp() + stats.endo->base.fn() : NAN);
    x["edge"]     = S0(o.edge);
    
    return x;
}

void VAlign::writeSummary(const FileName &file,
                          const FileName &gen,
                          const FileName &seq,
                          const VAlign::Stats &stats,
                          const VAlign::Options &o)
{
    const auto sumg2c = sum(stats.endo->covered);
    const auto sumg2l = sum(stats.endo->length);
    const auto sums2c = sum(stats.seqs->covered);
    const auto sums2l = sum(stats.seqs->length);
    
    A_ASSERT(sums2l >= sums2c);
    A_ASSERT(sumg2l >= sumg2c);
    
    const auto summary1 = "-------VarAlign Summary Statistics\n\n"
                          "       Reference annotation file: %1%\n"
                          "       Sample alignment file: %2%\n"
                          "       Sequin alignment file: %3%\n\n"
                          "       Edge: %24% bases per region\n\n"
                          "-------Alignments\n\n"
                          "       Sample: %4% (%5%%%)\n"
                          "       Sequin: %6% (%7%%%)\n"
                          "       Dilution:  %8$.4f%%\n\n"
                          "-------Reference regions\n\n"
                          "       Regions: %9% regions\n"
                          "       Regions: %10% bases\n\n"
                          "-------Alignment coverage of reference (Sample)\n\n"
                          "       *Nucleotide level\n"
                          "       Covered:     %20%\n"
                          "       Uncovered:   %21%\n"
                          "       Total:       %22%\n\n"
                          "       Sensitivity: %23$.4f\n\n"
                          "-------Alignment coverage of reference (Sequin)\n\n"
                          "       *Alignment level\n"
                          "       Inside regions:  %11%\n"
                          "       Outside regions: %12%\n\n"
                          "       Precision:       %13$.4f\n\n"
                          "       *Nucleotide level\n"
                          "       Covered:     %14%\n"
                          "       Uncovered:   %15%\n"
                          "       Erroneous:   %16%\n"
                          "       Total:       %17%\n\n"
                          "       Sensitivity: %18$.4f\n"
                          "       Precision:   %19$.4f\n";

    const auto summary2 = "-------VarAlign Summary Statistics\n\n"
                          "       Reference annotation file: %1%\n"
                          "       Sequin alignment file: %2%\n\n"
                          "       Edge: %15% bases per region\n\n"
                          "-------Alignments\n\n"
                          "       Synthetic: %3% \n\n"
                          "-------Reference regions\n\n"
                          "       Regions: %4% regions\n"
                          "       Regions: %5% bases\n\n"
                          "-------Comparison of alignments to annotation (Sequin)\n\n"
                          "       *Alignment level\n"
                          "       Inside regions:  %6%\n"
                          "       Outside regions: %7%\n\n"
                          "       Precision:      %8$.4f\n\n"
                          "       *Nucleotide level\n"
                          "       Covered:     %9%\n"
                          "       Uncovered:   %10%\n"
                          "       Erroneous:   %11%\n"
                          "       Total:       %12%\n\n"
                          "       Sensitivity: %13$.4f\n"
                          "       Precision:   %14$.4f\n";
    
    auto x = jsonD(gen, seq, stats, o);
    
    o.generate(file);
    o.writer->open(file);
    
    if (!gen.empty())
    {
        o.writer->write((boost::format(summary1) % x["ref"]      // 1
                                                 % x["endo"]     // 2
                                                 % x["seqs"]     // 3
                                                 % x["nEndo"]    // 4
                                                 % x["pEndo"]    // 5
                                                 % x["nSeqs"]    // 6
                                                 % x["pSeqs"]    // 7
                                                 % x["dilution"] // 8
                                                 % x["nRegs"]    // 9
                                                 % x["lRegs"]    // 10
                                                 % x["aTP"]      // 11
                                                 % x["aFP"]      // 12
                                                 % x["aPC"]      // 13
                                                 % x["bTP"]      // 14
                                                 % x["bFN"]      // 15
                                                 % x["bFP"]      // 16
                                                 % x["bTot"]     // 17
                                                 % x["bSN"]      // 18
                                                 % x["bPC"]      // 19
                                                 % x["eTP"]      // 20
                                                 % x["eFN"]      // 21
                                                 % x["eTot"]     // 22
                                                 % x["eSN"]      // 23
                                                 % x["edge"]     // 24
                         ).str());
    }
    else
    {
        o.writer->write((boost::format(summary2) % x["ref"]   // 1
                                                 % x["seqs"]  // 2
                                                 % x["nSeqs"] // 3
                                                 % x["nRegs"] // 4
                                                 % x["lRegs"] // 5
                                                 % x["aTP"]   // 6
                                                 % x["aFP"]   // 7
                                                 % x["aPC"]   // 8
                                                 % x["bTP"]   // 9
                                                 % x["bFN"]   // 10
                                                 % x["bFP"]   // 11
                                                 % x["bTot"]  // 12
                                                 % x["bSN"]   // 13
                                                 % x["bPC"]   // 14
                                                 % x["edge"]  // 15
                         ).str());
    }

    o.writer->close();
}

void VAlign::writeBQuins(const FileName &file,
                         const VAlign::Stats &stats,
                         const VAlign::Options &o)
{
#ifdef DEBUG_VALIGN
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Chrom" % "Position" % "Label").str());
    
    // For each chromosome...
    for (const auto &i : stats.seqs->data)
    {
        const auto &cID = i.first;
        
        // For each region in the chromosome...
        for (const auto &j : stats.seqs->inters.at(cID).data())
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
    const auto format = "%1%\t%2%\t%3%\t%4$.4f\t%5$.4f";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Length"
                                           % "Reads"
                                           % "Sn"
                                           % "Pc").str());

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

    for (const auto &i : stats.seqs->data)
    {
        for (const auto &j : i.second.afp)
        {
            o.writer->write((boost::format(format) % j % "FP").str());
        }
    }
    
    o.writer->close();
#endif
}

void VAlign::report(const FileName &endo, const FileName &seqs, const Options &o)
{
    const auto stats = analyze(endo, seqs, o);

    o.info("Generating statistics");
    
    /*
     * Generating VarAlign_summary.stats
     */
    
    writeSummary("VarAlign_summary.stats", endo, seqs, stats, o);

    /*
     * Generating VarAlign_sequins.tsv
     */
    
    writeQuins("VarAlign_sequins.tsv", stats, o);

    /*
     * Generating VarAlign_queries.stats (for debugging)
     */
    
    writeQueries("VarAlign_queries.stats", stats, o);

    /*
     * Generating VarAlign_rbase.stats (for debugging)
     */
    
    writeBQuins("VarAlign_rbase.stats", stats, o);
    
    /*
     * Generating VarAlign_stats.json
     */
    
    o.generate("VarAlign_stats.json");
    JSONWriter w(o.work);
    w.open("VarAlign_stats.json");
    w.write(jsonD(endo, seqs, stats, o));
    w.close();
}
