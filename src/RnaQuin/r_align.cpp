#include "RnaQuin/r_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;
using namespace std::placeholders;

// Defined in resources.cpp
extern FileName GTFRef();

static RAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;

    RAlign::Stats stats;

    stats.eInters = r.ueInters();
    stats.iInters = r.uiInters();

    assert(stats.eInters.size());
    assert(stats.iInters.size());
    assert(stats.eInters.size() == stats.iInters.size());

    for (const auto &i : stats.eInters)
    {
        stats.data[i.first];
    }

    assert(!stats.data.empty());
    return stats;
}

RAlign::Stats calculate(const RAlign::Options &o, std::function<void (RAlign::Stats &)> f)
{
    const auto &r = Standard::instance().r_trans;

    auto stats = init();
    
    // Parsing input files
    f(stats);

    o.info("Collecting statistics");
    
    for (const auto i : stats.eInters)
    {
        const auto &cID = i.first;
        const auto &x = stats.data.at(cID);

        const auto bs  = i.second.stats();
        const auto btp = bs.nonZeros;
        const auto bfp = x.bLvl.fp;
        const auto bfn = bs.length;
        
        if (Standard::isSynthetic(cID))
        {
            stats.sbm.tp() += btp;
            stats.sbm.fp() += bfp;
            stats.sbm.fn() += bfn;
            stats.sam.tp() += x.aLvl.m.tp();
            stats.sam.fp() += x.aLvl.m.fp();
            stats.sim.tp() += x.iLvl.tp();
            stats.sim.fp() += x.iLvl.fp();
            stats.sem.tp() += x.eLvl.tp();
        }
        else
        {
            stats.gbm.tp() += btp;
            stats.gbm.fp() += bfp;
            stats.gbm.fn() += bfn;
            stats.gam.tp() += x.aLvl.m.tp();
            stats.gam.fp() += x.aLvl.m.fp();
            stats.gim.tp() += x.iLvl.tp();
            stats.gim.fp() += x.iLvl.fp();
            stats.gem.tp() += x.eLvl.tp();
        }
    }
    
    stats.sem.fn() = r.countUExonSyn();
    stats.gem.fn() = r.countUExonGen();
    stats.sim.fn() = r.countUIntrSyn();
    stats.gim.fn() = r.countUIntrGen();

    return stats;
}

static void match(RAlign::Stats &stats, const ParserSAM::Data &align)
{
    Locus l;
    bool spliced;
    
    auto &x = stats.data.at(align.cID);
    
    if (align.spliced) { x.aLvl.spliced++; }
    else               { x.aLvl.normal++;  }
    
    // Is this a read a TP at the alignment level?
    auto isTP = true;
    
    while (isTP && align.nextCigar(l, spliced))
    {
        if (spliced)
        {
            // Can we find an exact match for the intron?
            const auto match = stats.iInters.at(align.cID).exact(l);
            
            if (match)
            {
                x.iLvl.tp()++;
            }
            else
            {
                isTP = false;
                x.iLvl.fp()++;
            }
        }
        else
        {
            // Can we find an overlapping match for the exon?
            const auto match = stats.eInters.at(align.cID).contains(l);
            
            if (match)
            {
                // This is the first alignment to the exon (we're only interested in the sensitivty)
                if (!match->stats().nonZeros)
                {
                    x.eLvl.tp()++;
                }
                
                const auto fp = match->map(l);
                assert(!fp);
            }
            else
            {
                isTP = false;
                
                // Can we find an overlapping match for the exon?
                const auto match = stats.eInters.at(align.cID).overlap(l);
                
                if (match)
                {
                    x.bLvl.fp += match->map(l);
                }
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

RAlign::Stats RAlign::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](RAlign::Stats &stats)
    {
        ParserSAM::parse(file, [&](const ParserSAM::Data &x, const ParserSAM::Info &info)
        {
            stats.update(x);

            if (!x.mapped)
            {
                return;
            }            
            else if (Standard::isSynthetic(x.cID) || Standard::isGenomic(x.cID))
            {
                match(stats, x);
            }
        });
    });
}

static Scripts summary()
{
    return "-------RnaAlign Summary Statistics\n\n"
           "       Input alignment file: %1%\n"
           "       Reference annotation file: %2%\n\n"
           "-------Number of alignments mapped to the synthetic chromosome and human genome\n\n"
           "       Synthetic: %3%\n"
           "       Genome:    %4%\n"
           "       Dilution:  %5$.2f\n\n"
           "-------Reference annotation (Synthetic)\n\n"
           "       Synthetic: %7% exons\n"
           "       Synthetic: %8% introns\n"
           "       Synthetic: %9% bases\n\n"
           "-------Reference annotation (Genome)\n\n"
           "       Genome: %10% exons\n"
           "       Genome: %11% introns\n"
           "       Genome: %12% bases\n\n"
           "-------Alignments (Synthetic)\n\n"
           "       Non-spliced: %13%\n"
           "       Spliced:     %14%\n"
           "       Total:       %15%\n\n"
           "-------Alignments (Genome)\n\n"
           "       Non-spliced: %16%\n"
           "       Spliced:     %17%\n"
           "       Total:       %18%\n\n"
           "-------Comparison of alignments to reference annotation (Synthetic)\n\n"
           "       *Exon level\n"
           "        Sensitivity: %19$.2f\n"
           "        Precision:   %20$.2f\n\n"
           "       *Intron level\n"
           "        Sensitivity: %21$.2f\n"
           "        Precision:   %22$.2f\n\n"
           "       *Base level\n"
           "        Sensitivity: %23$.2f\n"
           "        Precision:   %24$.2f\n\n"
           "-------Comparison of alignments to reference annotation (Genome)\n\n"
           "       *Exon level\n"
           "        Sensitivity: %28%\n"
           "        Precision:   %29%\n\n"
           "       *Intron level\n"
           "        Sensitivity: %30%\n"
           "        Precision:   %31%\n\n"
           "       *Base level\n"
           "        Sensitivity: %32%\n"
           "        Precision:   %33%\n\n";
}

static void generateSummary(const FileName &file,
                            const FileName &src,
                            const RAlign::Stats &stats,
                            const RAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto hasGeno = stats.data.size() > 1;
    
    #define CHECK(x) (hasGeno ? toString(x) : "-")
    
    o.writer->open(file);
    o.writer->write((boost::format(summary()) % src                      // 1
                                              % GTFRef()                 // 2
                                              % stats.n_syn              // 3
                                              % stats.n_gen              // 4
                                              % stats.dilution()         // 5
                                              % stats.n_unmap            // 6
                                              % r.countUExonSyn()        // 7
                                              % r.countUIntrSyn()        // 8
                                              % r.countLenSyn()          // 9
                                              % CHECK(r.countUExonGen()) // 10
                                              % CHECK(r.countUIntrGen()) // 11
                                              % CHECK(r.countLenGen())   // 12
                                              % ""//stats.countSpliceSyn()   // 13
                                              % ""//stats.countNormalSyn()   // 14
                                              % ""//(stats.countSpliceSyn() + stats.countNormalSyn())
                                              % ""//CHECK(stats.countSpliceGen())
                                              % ""//CHECK(stats.countNormalGen())
                                              % ""//CHECK(stats.countSpliceGen() + stats.countNormalGen())
                                              % ""//stats.s_esn        // 19
                                              % ""//stats.s_epc        // 20
                                              % ""//stats.s_isn        // 21
                                              % ""//stats.s_ipc        // 22
                                              % ""//stats.s_bsn        // 23
                                              % ""//stats.s_bpc        // 24
                                              % ""//stats.s_ems        // 25
                                              % ""//stats.s_ims        // 26
                                              % ""//stats.s_gms        // 27
                                              % ""//CHECK(stats.g_esn) // 28
                                              % ""//CHECK(stats.g_epc) // 29
                                              % ""//CHECK(stats.g_isn) // 30
                                              % ""//CHECK(stats.g_ipc) // 31
                                              % ""//CHECK(stats.g_bsn) // 32
                                              % ""//CHECK(stats.g_bpc) // 33
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file,
                       const FileName &src,
                       const RAlign::Stats &stats,
                       const RAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";

    const auto &data = stats.data.at(ChrT);
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Reads"
                                           % "Sn_exon"
                                           % "Sn_intron"
                                           % "Pc_intron"
                                           % "Sn_base"
                                           % "Pc_base").str());

//    for (const auto &i : data.overB.hist)
//    {
//        // Eg: R1_1
//        const auto &id = i.first;
//        
//        const auto &mb = data.geneB.at(id);
//        const auto &me = data.geneE.at(id);
//        const auto &mi = data.geneI.at(id);
//
//        const auto m = r.findGene(ChrT, id);
//        const auto reads = stats.s2r.at(id);
//        
//        assert(m);
//        
//        // Not all sequins have an intron...
//        if (mi.lNR)
//        {
//            o.writer->write((boost::format(format) % id
//                                                   % reads
//                                                   % me.sn()
//                                                   % me.pc()
//                                                   % mi.sn()
//                                                   % mi.pc()
//                                                   % mb.sn()
//                                                   % mb.pc()).str());
//        }
//        else
//        {
//            o.writer->write((boost::format(format) % id
//                                                   % reads
//                                                   % me.sn()
//                                                   % me.pc()
//                                                   % "--"
//                                                   % "--"
//                                                   % mb.sn()
//                                                   % mb.pc()).str());
//        }
//    }

    o.writer->close();
}

void RAlign::report(const FileName &file, const Options &o)
{
    const auto stats = RAlign::analyze(file, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating RnaAlign_summary.stats
     */
    
    o.analyze("RnaAlign_summary.stats");
    generateSummary("RnaAlign_summary.stats", file, stats, o);

    /*
     * Generating RnaAlign_quins.stats
     */
    
    o.analyze("RnaAlign_sequins.csv");
    writeQuins("RnaAlign_sequins.csv", file, stats, o);

    /*
     * Generating RnaAlign_report.pdf
     */
    
    o.report->open("RnaAlign_report.pdf");
    o.report->addTitle("RnaAlign");
    o.report->addFile("RnaAlign_summary.stats");
    o.report->addFile("RnaAlign_sequins.csv");
}