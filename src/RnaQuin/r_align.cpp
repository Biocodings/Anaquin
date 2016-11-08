#include "tools/errors.hpp"
#include "tools/gtf_data.hpp"
#include "RnaQuin/r_align.hpp"
#include "RnaQuin/RnaQuin.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern FileName GTFRef();

#ifdef ANAQUIN_DEBUG
static std::ofstream __iWriter__;
static std::ofstream __bWriter__;
static std::ofstream __rWriter__;
#endif

static void writeIntron(const ChrID &cID, const Locus &l, const GeneID &gID, const Label &label)
{
#ifdef ANAQUIN_DEBUG
    __iWriter__ << cID << "\t" << l.start << "-" << l.end << "\t" << gID << "\t" << label << "\n";
#endif
}

static void writeBase(const ChrID &cID, const Locus &l, const Label &label)
{
#ifdef ANAQUIN_DEBUG
    __bWriter__ << cID << "\t" << l.start << "-" << l.end << "\t" << label << "\n";
#endif
}

static RAlign::Stats init()
{
    const auto &r = Standard::instance().r_rna;

    RAlign::Stats stats;

    stats.iInters = r.uiInters();

    /*
     * It's important to use meInters() rather than ueInters(). Due to alternative splicing, two
     * different transcripts can give overlapping exons. We don't know where exactly the read
     * come from. The exon in transcript 1? The exon in transcript 2? We don't have the information.
     * But we can construct non-overlapping (merged) exon regions and use that to calculate statistics
     * such as base-level sensitivity.
     *
     * It's important to note exon sensitivity is not possible here (again, due to alternative splicing).
     * The exon intervals are actually merged intervals.
     */

    stats.eInters = r.meInters(Strand::Either);

    A_CHECK(stats.eInters.size(), "stats.eInters.size()");
    A_CHECK(stats.iInters.size(), "stats.iInters.size()");
    A_CHECK(stats.eInters.size() == stats.iInters.size(), "stats.eInters.size() == stats.iInters.size()");
    
    for (const auto &i : stats.eInters)
    {
        const auto &cID = i.first;
        
        // Number of unique reference exons
        stats.data[cID].eLvl.nr() = r.countUExon(cID);
        
        // Number of unique reference introns
        stats.data[cID].iLvl.m.nr() = r.countUIntr(cID);

        /*
         * We'd like to know the length of the chromosome but we don't have the information.
         * It doesn't matter because we can simply initalize it to the the maximum possible.
         * The data structure in MergedInterval will be efficient not to waste memory.
         */

        MergedInterval *mi = new MergedInterval(cID, Locus(1, std::numeric_limits<Base>::max()));
        stats.data[cID].bLvl.fp = std::shared_ptr<MergedInterval>(mi);
        
        A_CHECK(stats.data[cID].eLvl.nr(), "stats.data[cID].eLvl.nr()");
    }

    A_CHECK(!stats.data.empty(), "!stats.data.empty()");
    return stats;
}

RAlign::Stats calculate(const RAlign::Options &o, std::function<void (RAlign::Stats &)> f)
{
    const auto &r = Standard::instance().r_rna;

    auto stats = init();

#ifdef ANAQUIN_DEBUG
    __bWriter__.open(o.work + "/RnaAlign_qbase.txt");
    __iWriter__.open(o.work + "/RnaAlign_qintrs.txt");
    __rWriter__.open(o.work + "/RnaAlign_reads.txt");
#endif

    // Parsing input files
    f(stats);

#ifdef ANAQUIN_DEBUG
    __iWriter__.close();
    __bWriter__.close();
#endif

    o.info("Collecting statistics");
    
    o.logInfo("Reference chromsomes: == " + std::to_string(stats.data.size()));
    o.logInfo("Exon intervals: " + std::to_string(stats.eInters.size()));
    
    // For each reference chromosome...
    for (const auto i : stats.eInters)
    {
        const auto &cID = i.first;
        const auto &x = stats.data.at(cID);

        const auto bs = i.second.stats();
        
        /*
         * Calculating statistics for alignments
         */

        // Number of alignments
        const auto atp = x.aLvl.m.tp();
        
        // Number of alignments
        const auto afp = x.aLvl.m.fp();

        /*
         * Calculating statistics for unique exons
         */
        
        const auto es = stats.eInters.at(cID).stats();
        
        /*
         * Calculating statistics for bases
         */

        const auto bfp = x.bLvl.fp->stats().nonZeros;
        const auto btp = es.nonZeros;
        const auto bfn = bs.length - bs.nonZeros;

        /*
         * Calculating statistics for unique introns
         */

        const auto is = stats.iInters.at(cID).stats();
        
        // Number of unique introns exactly detected
        const auto itp = is.f;
        
        // Number of unique introns not detected
        const auto ifn = is.n - is.f;

        // Number of unique introns predicted but not detected in the reference
        const auto ifp = x.iLvl.fp.size();
        
        /*
         * Aggregating statistics for synthetic chromosomes and genomic chromosomes
         */
        
        if (isRnaQuin(cID))
        {
            stats.sbm.tp() += btp;
            stats.sbm.fp() += bfp;
            stats.sbm.fn() += bfn;

            stats.sam.tp() += atp;
            stats.sam.fp() += afp;

            stats.sem.tp() += x.eLvl.tp();
            stats.sem.fn() += (r.countUExon(cID) - x.eLvl.tp());
            
            stats.sim.fp() += ifp;
            stats.sim.tp() += itp;
            stats.sim.fn() += ifn;

            stats.sn += x.aLvl.normal;
            stats.ss += x.aLvl.spliced;
        }
        else
        {
            stats.gbm.tp() += btp;
            stats.gbm.fp() += bfp;
            stats.gbm.fn() += bfn;

            stats.gam.tp() += atp;
            stats.gam.fp() += afp;

            stats.gem.tp() += x.eLvl.tp();
            stats.gem.fn() += (r.countUExon(cID) - x.eLvl.tp());

            stats.gim.fp() += ifp;
            stats.gim.tp() += itp;
            stats.gim.fn() += ifn;
            
            stats.gn += x.aLvl.normal;
            stats.gs += x.aLvl.spliced;
        }
    }
    
    stats.sem.fn() = r.countUExonSyn();
    stats.gem.fn() = r.countUExonGen();

    return stats;
}

static void match(RAlign::Stats &stats, const ParserSAM::Info &info, ParserSAM::Data &align)
{
    static Locus l;
    static bool spliced;

    auto &x = stats.data.at(align.cID);

    if (info.skip)
    {
        x.aLvl.spliced++;
    }
    else
    {
        // Indels (ins+del) are also counted as "normal"
        x.aLvl.normal++;
    }

    // This'll be set to false whenever there is a mismatch
    bool isTP = true;
    
    GeneID gID = "";

    // Check all cigar blocks...
    while (align.nextCigar(l, spliced))
    {
        if (spliced)
        {
            // Can we find an exact match for the intron?
            auto match = stats.iInters.at(align.cID).exact(l);
            
            if (match)
            {
                // We'll use it to calculate sensitivty at the intron level
                match->map(l);

                writeIntron(align.cID, l, match->gID(), "TP");
            }
            else
            {
                x.iLvl.fp.insert(l);
                isTP = false;

                writeIntron(align.cID, l, "", "FP");
            }
        }
        else
        {
            // Can we find an contained match for the exon?
            const auto match = stats.eInters.at(align.cID).contains(l);
            
#ifdef DEBUG_ANAQUIN
            if (ms.size() > 1)
            {
                o.logInfo("Mult1: " + align.cID + "\t" + match->tID() + "\t" + std::to_string(l.start) + "-" std::to_string(l.end));
                o.logInfo("Mult2: " + align.cID + "\t" + match->tID() + "\t" + std::to_string(match->l.start) + "-" std::to_string(match->l.end));
            }
#endif
            
            if (match)
            {
                // We'll need it for calculating sensitivity at the base level
                match->map(l);
                
                gID = match->gID();

                writeBase(align.cID, l, "TP");
            }
            else
            {
                // Can we find an overlapping match for the exon?
                const auto match = stats.eInters.at(align.cID).overlap(l);

                if (match)
                {
                    match->map(l);
                    
                    // Gap to the left?
                    if (l.start < match->l().start)
                    {
                        const auto gap = Locus(l.start, match->l().start-1);
                        
                        x.bLvl.fp->map(gap);
                        
                        writeBase(align.cID, gap, "FP");
                    }
                    
                    // Gap to the right?
                    if (l.end > match->l().end)
                    {
                        const auto gap = Locus(match->l().end+1, l.end);
                        
                        x.bLvl.fp->map(gap);
                        
                        writeBase(align.cID, gap, "FP");
                    }
                }
                else
                {
                    // The entire locus is outside of the reference region
                    x.bLvl.fp->map(l);
                    
                    writeBase(align.cID, l, "FPO");
                }
                
                // The alignment is overlapping, thus it's a FP
                isTP = false;
            }
        }
    }
    
    if (isTP)
    {
        x.aLvl.m.tp()++;

        A_CHECK(!gID.empty(), "!gID.empty()");
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
            if (info.p.i && !(info.p.i % 1000000))
            {
                o.wait(std::to_string(info.p.i));
            }

            // Don't count for multiple alignments
            if (!x.mapped || x.isPrimary)
            {
#ifdef ANAQUIN_DEBUG
                if (x.mapped && x.cID != ChrIS)
                    __rWriter__ << x.name << "\n";
#endif
                stats.update(x, isRnaQuin);
            }

            if (!x.mapped)
            {
                return;
            }
            else if (isRnaQuin(x.cID) || Standard::isGenomic(x.cID))
            {
                match(stats, info, x);
            }
            else
            {
                o.logWarn("Ignore: " + x.name + "  " + x.cID);
            }
        });
    });
}

static Scripts summary()
{
    return "-------RnaAlign Summary Statistics\n\n"
           "       Input alignment file: %1%\n"
           "       Reference annotation file: %2%\n\n"
           "-------Number of alignments mapped to the synthetic chromosome and genome\n\n"
           "       Synthetic: %3%\n"
           "       Genome:    %4%\n"
           "       Dilution:  %5$.3f\n\n"
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
           "       Spliced:     %14%\n\n"
           "-------Alignments (Genome)\n\n"
           "       Non-spliced: %15%\n"
           "       Spliced:     %16%\n\n"
           "-------Comparison of alignments to reference annotation (Synthetic)\n\n"
           "       *Intron level\n"
           "        Sensitivity: %17$.2f\n"
           "        Precision:   %18$.2f\n\n"
           "       *Base level\n"
           "        Sensitivity: %19$.2f\n"
           "        Precision:   %20$.2f\n\n"
           "-------Comparison of alignments to reference annotation (Genome)\n\n"
           "       *Intron level\n"
           "        Sensitivity: %21$.2f\n"
           "        Precision:   %22$.2f\n\n"
           "       *Base level\n"
           "        Sensitivity: %23$.2f\n"
           "        Precision:   %24$.2f\n";
}

static void generateSummary(const FileName &file,
                            const FileName &src,
                            const RAlign::Stats &stats,
                            const RAlign::Options &o)
{
    const auto &r = Standard::instance().r_rna;
    const auto hasGeno = stats.data.size() > 1;
    
    #define G(x) (hasGeno ? toString(x) : "-")
    
    o.writer->open(file);
    o.writer->write((boost::format(summary()) % src                  // 1
                                              % GTFRef()             // 2
                                              % stats.nSyn           // 3
                                              % stats.nGen           // 4
                                              % stats.dilution()     // 5
                                              % stats.nNA            // 6
                                              % r.countUExonSyn()    // 7
                                              % r.countUIntrSyn()    // 8
                                              % r.countLenSyn()      // 9
                                              % G(r.countUExonGen()) // 10
                                              % G(r.countUIntrGen()) // 11
                                              % G(r.countLenGen())   // 12
                                              % stats.sn             // 13
                                              % stats.ss             // 14
                                              % G(stats.gn)          // 15
                                              % G(stats.gs)          // 16
                                              % stats.sim.sn()       // 17
                                              % stats.sim.pc()       // 18
                                              % stats.sbm.sn()       // 19
                                              % stats.sbm.pc()       // 20
                                              % stats.gim.sn()       // 21
                                              % stats.gim.pc()       // 22
                                              % stats.gbm.sn()       // 23
                                              % stats.gbm.pc()       // 24
                     ).str());
    o.writer->close();
}

static void writeBQuins(const FileName &file,
                        const FileName &src,
                        const RAlign::Stats &stats,
                        const RAlign::Options &o)
{
#ifdef ANAQUIN_DEBUG
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ChrID" % "Position" % "Label").str());
    
    for (const auto &i : stats.data)
    {
        const auto &cID = i.first;
        
        for (const auto &j : stats.eInters.at(cID).data())
        {
            for (const auto &k : j.second._data)
            {
                const auto pos = (toString(k.second.start) + "-" + toString(k.second.end));
                o.writer->write((boost::format(format) % cID % pos % "TP").str());
            }
        }
    }
    
    o.writer->close();
#endif
}

static void writeIQuins(const FileName &file,
                        const FileName &src,
                        const RAlign::Stats &stats,
                        const RAlign::Options &o)
{
#ifdef ANAQUIN_DEBUG
    const auto format = "%1%\t%2%\t%3%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "ChrID" % "Position" % "Label").str());
    
    for (const auto &i : stats.data)
    {
        const auto &cID = i.first;
        
        for (const auto &j : stats.iInters.at(cID).data())
        {
            auto is = j.second.stats();
            
            A_CHECK(is.nonZeros == 0 || is.nonZeros == is.length, "is.nonZeros == 0 || is.nonZeros == is.length");
            
            const auto pos = (toString(j.second.l().start) + "-" + toString(j.second.l().end));
            
            if (is.nonZeros)
            {
                o.writer->write((boost::format(format) % cID % pos % "TP").str());
            }
            else
            {
                o.writer->write((boost::format(format) % cID % pos % "FN").str());
            }
        }
    }
    
    o.writer->close();
#endif
}

static void writeQuins(const FileName &file,
                       const FileName &src,
                       const RAlign::Stats &stats,
                       const RAlign::Options &o)
{
    const auto &r = Standard::instance().r_rna;

    const auto h2g = r.histGene();
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";

    o.writer->open(file);
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Reads"
                                           % "SnIntron"
                                           % "SnBase").str());

    for (const auto &i : stats.data)
    {
        const auto &cID = i.first;
        
        if (isRnaQuin(cID))
        {
            std::map<GeneID, Confusion> bm, im;

#ifdef ANAQUIN_DEBUG
            std::map<GeneID, Confusion> em;
            
            /*
             * Calculating exon statistics for the genes
             */
            
            for (const auto &j : stats.eInters.at(cID).data())
            {
                const auto &gID = j.second.gID();
                
                // Statistics for the exon within the gene
                const auto is = j.second.stats();
                
                if (!is.nonZeros)
                {
                    o.logInfo("Exon: (FN): " + gID + " " + std::to_string(j.second.l().start) + "-" + std::to_string(j.second.l().end));
                    em[gID].fn()++;
                }
                else
                {
                    em[gID].tp()++;
                }
            }
#endif
            
            /*
             * Calculating intron statistics for the genes
             */
            
            for (const auto &j : stats.iInters.at(cID).data())
            {
                const auto &gID = j.second.gID();

                // Statistics for the intron within the gene
                const auto is = j.second.stats();
                
                A_CHECK(is.nonZeros == 0 || is.nonZeros == is.length, "is.nonZeros == 0 || is.nonZeros == is.length");
                
                if (!is.nonZeros)
                {
                    o.logInfo("Intron: (FN): " + gID + " " + std::to_string(j.second.l().start) + "-" + std::to_string(j.second.l().end));
                    im[gID].fn()++;
                }
                else
                {
                    im[gID].tp()++;
                }
            }
            
            /*
             * Calculating base statistics for the genes
             */
            
            for (const auto &j : stats.eInters.at(cID).data())
            {
                const auto &gID = j.second.gID();
                
                // Statistics for the bases within the gene
                const auto bs = j.second.stats();

                bm[gID].tp() += bs.nonZeros;
                bm[gID].fn() += bs.length - bs.nonZeros;
            }
            
            // For every gene in the reference
            for (const auto &j : h2g.at(cID))
            {
                const auto &data = i.second;
                
                // Eg: R1_1
                const auto &gID = j.first;
                
                // Number of reads aligned
                const auto reads = data.g2r.count(gID) ? data.g2r.at(gID) : 0;

                // Sensitivity at the intron level
                const auto isn = im.count(gID) ? std::to_string(im.at(gID).sn()) : "-";
                
                o.writer->write((boost::format(format) % gID
                                                       % r.findGene(cID, gID)->l.length()
                                                       % reads
                                                       % isn
                                                       % bm.at(gID).sn()).str());
            }
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
    
    o.generate("RnaAlign_summary.stats");
    generateSummary("RnaAlign_summary.stats", file, stats, o);

    /*
     * Generating RnaAlign_sequins.csv
     */
    
    o.generate("RnaAlign_sequins.csv");
    writeQuins("RnaAlign_sequins.csv", file, stats, o);

    /*
     * Generating RnaAlign_rintrs.txt
     */
    
    writeIQuins("RnaAlign_rintrs.txt", file, stats, o);

    /*
     * Generating RnaAlign_rbase.txt
     */
    
    writeBQuins("RnaAlign_rbase.txt", file, stats, o);
}
