#include "VarQuin/v_align.hpp"
#include "VarQuin/VarQuin.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

static VAlign::Stats init()
{
    const auto &r = Standard::instance().r_var;

    VAlign::Stats stats;
    
    stats.hist   = r.hist();
    stats.inters = r.inters();
    
    assert(!stats.hist.empty());
    assert(!stats.inters.empty());

    return stats;
}

static void classifyAlign(VAlign::Stats &stats, const Alignment &align)
{
    Base lGaps, rGaps;
    
    auto f = [&](Interval *m)
    {
        m->map(align.l, &lGaps, &rGaps);

        const auto covered = (align.l.length() - lGaps - rGaps);
        
        stats.data[align.cID].lGaps[m->name()] += lGaps;
        stats.data[align.cID].lGaps[m->name()] += rGaps;
        stats.data[align.cID].align[m->name()] += covered;
        
        assert(covered >= 0);
        assert(align.l.length() > lGaps);
        assert(align.l.length() > rGaps);
    };
    
    // Does the read aligned within a gene (or a region)?
    const auto m = stats.inters.at(align.cID).contains(align.l);

    if (m)
    {
        f(m);
        assert(lGaps == 0 && rGaps == 0);

        stats.data[align.cID].tp++;
        //stats.data[align.cID].gtp[m->name()]++;
        stats.hist.at(align.cID).at(m->name())++;
    }
    else
    {
        // Can we at least match by overlapping?
        const auto m = stats.inters[align.cID].overlap(align.l);

        if (m)
        {
            f(m);
            assert(lGaps != 0 || rGaps != 0);

            stats.data[align.cID].fp++;
            //stats.data[align.cID].gfp[m->name()]++;
            stats.data[align.cID].afp.push_back(align.name);
        }
        else if (Standard::isSynthetic(align.cID))
        {
            stats.data[align.cID].fp++;

            /*
             * The read is not aligned within the reference regions. We don't know whether this is
             * a TP or FP.
             */
        }
    }
}

VAlign::Stats VAlign::analyze(const FileName &file, const Options &o)
{
    auto stats = init();
    o.analyze(file);

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &info)
    {
        if (!align.i && info.p.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        if (align.spliced)
        {
            o.warn("Spliced alignment detected");
        }

        stats.update(align);
        
        if (!align.mapped)
        {
            return;
        }
        else if (Standard::isSynthetic(align.cID))
        {
            classifyAlign(stats, align);
        }
        else if (Standard::isGenomic(align.cID))
        {
            classifyAlign(stats, align);
        }
    });

    /*
     * -------------------- Calculating statistics --------------------
     */

    Base stp = 0;
    Base sfp = 0;
    Base gtp = 0;
    Base gfp = 0;

    for (const auto &i : stats.hist)
    {
        for (const auto &j : i.second)
        {
            const auto &cID = i.first;
            const auto &gID = j.first;
         
            // Reads aligned to the gene
            stats.g2r[gID] = j.second;

            const auto isSyn = Standard::isSynthetic(cID);
            
            const auto m = stats.inters.at(cID).find(gID);
            assert(m);
            
            // Statistics for the gene (created by the interval)
            const auto ms = m->stats();

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
                stats.g2c[gID] = ms.covered();

                // Sensitivty for the gene
                stats.g2s[gID] = static_cast<Proportion>(stats.g2c.at(gID)) / stats.g2l.at(gID);
            }
            
            assert(stats.s2c[gID] <= stats.s2l[gID]);
            assert(stats.g2c[gID] <= stats.g2l[gID]);

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
            
            assert(stp >= 0);
            assert(sfp >= 0);

            stats.g2p[gID] = bpc;
        }
    }

    assert(!stats.g2r.empty());
    assert(!stats.g2s.empty());
    assert(!stats.s2l.empty());
    assert(!stats.s2c.empty());

    assert(stats.s2l.size() == stats.s2c.size());
    assert(stats.g2r.size() == stats.g2s.size());

    stats.spc = static_cast<Proportion>(stp) / (stp + sfp);
    stats.ssn = static_cast<Proportion>(sum(stats.s2c)) / sum(stats.s2l);

    stats.gpc = static_cast<Proportion>(gtp) / (gtp + gfp);
    stats.gsn = static_cast<Proportion>(sum(stats.g2c)) / sum(stats.g2l);

    assert(stats.spc >= 0.0 && stats.spc <= 1.0);
    assert(isnan(stats.gpc) || (stats.gpc >= 1.0 && stats.gpc <= 1.0));

    assert(stats.ssn >= 0.0 && stats.ssn <= 1.0);
    assert(isnan(stats.gsn) || (stats.gsn >= 1.0 && stats.gsn <= 1.0));
    
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
    assert(sumg2c >= sumg2l);

    const auto summary = "-------VarAlign Summary Statistics\n\n"
                         "       Reference annotation file: %1%\n"
                         "       User alignment file: %2%\n\n"
                         "-------Alignments\n\n"
                         "       Unmapped:  %3%\n"
                         "       Synthetic: %4% (%5%)\n"
                         "       Genome:    %6% (%7%)\n"
                         "       Dilution:  %8%\n\n"    
                         "-------Reference annotation (Synthetic)\n\n"
                         "       Synthetic: %9% genes\n"
                         "       Synthetic: %10% bases\n\n"
                         "-------Reference annotation (Genome)\n\n"
                         "       Genome: %11% genes\n"
                         "       Genome: %12% bases\n\n"
                         "-------Comparison of alignments to annotation (Synthetic)\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %13$.2f\n"
                         "       Uncovered:   %14$.2f\n"
                         "       Total:       %15$.2f\n"
                         "       Sensitivity: %16$.2f\n"
                         "       Precision:   %17$.2f\n\n"
                         "-------Comparison of alignments to annotation (Genome)\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %18$.2f\n"
                         "       Uncovered:   %19$.2f\n"
                         "       Total:       %20$.2f\n"
                         "       Sensitivity: %21$.2f\n"
                         "       Precision:   %22$.2f";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % BedRef()                // 1
                                            % src                     // 2
                                            % stats.n_unmap           // 3
                                            % stats.n_syn             // 4
                                            % (100 * stats.synProp()) // 5
                                            % stats.n_gen             // 6
                                            % (100 * stats.genProp()) // 7
                                            % stats.dilution()        // 8
                                            % r.countGeneSyn()        // 9
                                            % r.countBaseSyn()        // 10
                                            % r.countGeneGen()        // 11
                                            % r.countBaseGen()        // 12
                                            % sums2c                  // 13
                                            % (sums2l - sums2c)       // 14
                                            % sums2l                  // 15
                                            % stats.ssn               // 16
                                            % stats.spc               // 17
                                            % sumg2c                  // 18
                                            % (sumg2l - sumg2c)       // 19
                                            % sumg2l                  // 20
                                            % stats.gsn               // 21
                                            % stats.gpc               // 22
                     ).str());
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

    for (const auto &i : stats.hist)
    {
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

    for (const auto &j : stats.data.at(ChrT).afp)
    {
        o.writer->write((boost::format(format) % j % "FP").str());
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
     * Generating VarAlign_report.pdf
     */

    o.report->open("VarAlign_report.pdf");
    o.report->addTitle("VarAlign");
    o.report->addFile("VarAlign_summary.stats");
    o.report->addFile("VarAlign_sequins.csv");
    o.report->addFile("VarAlign_queries.stats");
}