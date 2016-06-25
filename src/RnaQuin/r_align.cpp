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
        const auto &cID = i.first;
        
        stats.data[cID];
        stats.data[cID].eLvl.nr() = r.countUExon(cID);
        stats.data[cID].iLvl.sn.nr() = r.countUIntr(cID);
        
        assert(stats.data[cID].eLvl.nr());
        assert(stats.data[cID].eLvl.nr() == 1 || stats.data[cID].iLvl.sn.nr());
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
        //const auto bfp = x.bLvl.fp();
        const auto bfn = bs.length - bs.nonZeros;
        
        if (Standard::isSynthetic(cID))
        {
            stats.sbm.tp() += btp;
            //stats.sbm.fp() += bfp;
            stats.sbm.fn() += bfn;
            stats.sam.tp() += x.aLvl.m.tp();
            stats.sam.fp() += x.aLvl.m.fp();

            stats.sem.tp() += x.eLvl.tp();
            stats.sem.fn() += (r.countUExon(cID) - x.eLvl.tp());

            stats.simpc += x.iLvl.pc;
            stats.simsn.tp() += x.iLvl.sn.tp();
            stats.simsn.fn() += (r.countUIntr(cID) - x.iLvl.sn.tp());

            stats.sn += x.aLvl.normal;
            stats.ss += x.aLvl.spliced;
        }
        else
        {
            stats.gbm.tp() += btp;
            //stats.gbm.fp() += bfp;
            stats.gbm.fn() += bfn;
            stats.gam.tp() += x.aLvl.m.tp();
            stats.gam.fp() += x.aLvl.m.fp();

            stats.gem.tp() += x.eLvl.tp();
            stats.gem.fn() += (r.countUExon(cID) - x.eLvl.tp());

            stats.gimpc += x.iLvl.pc;
            stats.gimsn.tp() += x.iLvl.sn.tp();
            stats.gimsn.fn() += (r.countUIntr(cID) - x.iLvl.sn.tp());
            
            stats.gn += x.aLvl.normal;
            stats.gs += x.aLvl.spliced;
        }
    }
    
    stats.sem.fn() = r.countUExonSyn();
    stats.gem.fn() = r.countUExonGen();

    return stats;
}

static void match(RAlign::Stats &stats, ParserSAM::Data &align)
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
                // Is this is the first alignment to the intron?
                if (!match->stats().nonZeros)
                {
                    x.iLvl.pc.tp()++;
                    x.iLvl.sn.tp()++;
                }
                
                match->map(l);
            }
            else
            {
                isTP = false;
                x.iLvl.pc.fp()++;
            }
        }
        else
        {
            // Can we find an overlapping match for the exon?
            const auto match = stats.eInters.at(align.cID).contains(l);
            
            if (match)
            {
                // Is this the first alignment to the exon?
                if (!match->stats().nonZeros)
                {
                    x.eLvl.tp()++;
                }
                
                match->map(l);
                //x.bLvl.tp() += l.length();
            }
            else
            {
                isTP = false;
                
                // Can we find an overlapping match for the exon?
                const auto match = stats.eInters.at(align.cID).overlap(l);
                
                if (match)
                {
                    match->map(l);
                }
            }
        }
    }

    if (isTP)
    {
        x.aLvl.m.tp()++;
        
        const auto &r = Standard::instance().r_trans;
        
        // By definition, our last block must be an alignment to an exon
        const auto &gID = r.findGene(align.cID, l, MatchRule::Overlap)->id();
        
        x.g2r[gID]++;
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
        ParserSAM::parse(file, [&](ParserSAM::Data &x, const ParserSAM::Info &info)
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
           "        Sensitivity: %19$.2f\n\n"
           "       *Intron level\n"
           "        Sensitivity: %20$.2f\n"
           "        Precision:   %21$.2f\n\n"
           "       *Base level\n"
           "        Sensitivity: %22$.2f\n"
           "        Precision:   %23$.2f\n\n"
           "-------Comparison of alignments to reference annotation (Genome)\n\n"
           "       *Exon level\n"
           "        Sensitivity: %24$.2f\n\n"
           "       *Intron level\n"
           "        Sensitivity: %25$.2f\n"
           "        Precision:   %26$.2f\n\n"
           "       *Base level\n"
           "        Sensitivity: %27$.2f\n"
           "        Precision:   %28$.2f\n\n";
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
                                              % stats.sn                 // 13
                                              % stats.ss                 // 14
                                              % (stats.sn + stats.ss)    // 15
                                              % stats.gn                 // 16
                                              % stats.gs                 // 17
                                              % (stats.gn + stats.gs)    // 18
                                              % stats.sem.sn()           // 19
                                              % stats.simsn.sn()         // 20
                                              % stats.simpc.pc()         // 21
                                              % stats.sbm.sn()           // 22
                                              % ""                       // 23
                                              % stats.gem.sn()           // 24
                                              % stats.gimsn.sn()         // 25
                                              % stats.gimpc.pc()         // 26
                                              % stats.gbm.sn()           // 27
                                              % ""                       // 28
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file,
                       const FileName &src,
                       const RAlign::Stats &stats,
                       const RAlign::Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";

    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Reads"
                                           % "Sn_exon"
                                           % "Sn_intron"
                                           % "Pc_intron"
                                           % "Sn_base").str());

    for (const auto &i : stats.data)
    {
        if (Standard::isSynthetic(i.first))
        {
            const auto &data = i.second;
            
            // Eg: R1_1
            const auto &id = i.first;
            
            // Number of reads aligned
            const auto reads = data.g2r.count(id) ? data.g2r.at(id) : 0;
            
            stats.eInters.at(i.first).find(id)->id();
            
            o.writer->write((boost::format(format) % id
                                                   % reads
                                                   % data.eLvl.sn()
                                                   % data.iLvl.sn.sn()
                                                   % data.iLvl.pc.pc()
                                                   % data.bLvl.sn()).str());
        }
    }

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