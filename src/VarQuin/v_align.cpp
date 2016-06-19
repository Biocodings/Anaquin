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

        stats.data[align.cID].lGaps[m->name()] += lGaps;
        stats.data[align.cID].lGaps[m->name()] += rGaps;
        stats.data[align.cID].align[m->name()] += (align.l.length() - lGaps - rGaps);
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
        else
        {
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

    auto stp = 0;
    auto sfp = 0;
    auto gtp = 0;
    auto gfp = 0;

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
            
            // Statistics for the gene
            const auto ms = m->stats();
            
            if (isSyn)
            {
                stats.s2l[gID] = ms.length;
                stats.s2c[gID] = ms.covered();

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

            // TP at the base level
            const auto btp = stats.data.at(i.first).align.at(gID);
            
            // FP at the base level
            const auto bfp = (stats.data.at(i.first).lGaps.count(gID) ? stats.data.at(i.first).lGaps.at(gID) : 0)
                                    +
                             (stats.data.at(i.first).rGaps.count(gID) ? stats.data.at(i.first).rGaps.at(gID) : 0);
            
            // Precison at the base level
            const auto bpc = static_cast<Proportion>(btp) / (btp + bfp);
            
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
    }

    assert(!stats.g2r.empty());
    assert(!stats.g2s.empty());
    assert(!stats.s2l.empty());
    assert(!stats.s2c.empty());

    assert(stats.s2l.size() == stats.s2c.size());
    assert(stats.g2r.size() == stats.g2s.size());

    stats.spc = static_cast<Proportion>(stp) / (stp + sfp);
    stats.gpc = static_cast<Proportion>(gtp) / (gtp + gfp);
    stats.ssn = static_cast<Proportion>(sum(stats.s2c)) / sum(stats.s2l);
    stats.gsn = static_cast<Proportion>(sum(stats.g2c)) / sum(stats.g2l);

    return stats;
}

static void writeSummary(const FileName &file, const FileName &src, const VAlign::Stats &stats, const VAlign::Options &o)
{
    extern FileName BedRef();

    const auto &r = Standard::instance().r_var;

    const auto sums2c = sum(stats.s2c);
    const auto sums2l = sum(stats.s2l);

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
                         "       *Region level\n"
                         "       Covered:     %13%\n"
                         "       Uncovered:   %14%\n"
                         "       Total:       %15%\n"
                         "       Sensitivity: %16%\n"
                         "       Precision:   %17%\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %18%\n"
                         "       Uncovered:   %19%\n"
                         "       Total:       %20%\n"
                         "       Sensitivity: %21%\n"
                         "       Precision:   %22%\n\n"
                         "-------Comparison of alignments to annotation (Genome)\n\n"
                         "       *Region level\n"
                         "       Covered:     %23%\n"
                         "       Uncovered:   %24%\n"
                         "       Total:       %25%\n"
                         "       Sensitivity: %26%\n"
                         "       Precision:   %27%\n\n"
                         "       *Nucleotide level\n"
                         "       Covered:     %28%\n"
                         "       Uncovered:   %29%\n"
                         "       Total:       %30%\n"
                         "       Sensitivity: %31%\n"
                         "       Precision:   %32%";

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
                                            % "????"                  // 13
                                            % "????"                  // 14
                                            % "????"                  // 15
                                            % "????"                  // 16
                                            % "????"                  // 17
                                            % sums2c                  // 18
                                            % (sums2l - sums2c)       // 19
                                            % sums2l                  // 20
                                            % stats.ssn               // 21
                                            % stats.spc               // 22
                                            % ""                      // 23
                                            % ""                      // 24
                                            % ""                      // 25
                                            % ""                      // 26
                                            % ""                      // 27
                                            % ""                      // 28
                                            % ""                      // 29
                                            % ""                      // 30
                                            % stats.gsn               // 31
                                            % stats.gpc               // 32
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

    for (const auto &i : stats.hist.at(ChrT))
    {
        o.writer->write((boost::format(format) % i.first
                                               % stats.s2l.at(i.first)
                                               % stats.g2r.at(i.first)
                                               % stats.g2s.at(i.first)
                                               % stats.g2p.at(i.first)).str());
    }

    o.writer->close();
}

static void writeQueries(const FileName &file, const VAlign::Stats &stats, const VAlign::Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "reads" % "label").str());

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