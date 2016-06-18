#include "RnaQuin/r_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;
using namespace std::placeholders;

// Defined in resources.cpp
extern Scripts PlotScatter();

// Defined in resources.cpp
extern FileName GTFRef();

// Internal implementation
typedef std::function<void (RAlign::Stats &)> Functor;

typedef RAlign::Stats::AlignMetrics   Metrics;
typedef RAlign::Stats::MissingMetrics MissingMetrics;

template <typename T> void initT(const ChrID &cID, T &t)
{
    const auto &r = Standard::instance().r_trans;

    // Initalize the distributions
    t.overB.hist = t.histE = t.histI = r.histGene()[cID];
    
    assert(!t.histE.empty());
    assert(!t.histI.empty());
    assert(!t.overB.hist.empty());
    
    /*
     * Initialize intervals for exons and introns
     */
    
    t.eInters = r.exonInters(cID);
    t.iInters = r.intronInters(cID);
    
    assert(t.eInters.size());
    assert(t.iInters.size());
    
    /*
     * Initialize overall statistics
     */
    
    for (const auto &i : t.histE)
    {
        t.geneB[i.first];
        t.geneE[i.first];
        t.geneI[i.first];
    }

    /*
     * Initialize statistics for reference exons
     */
    
    for (const auto &i : t.eInters.data())
    {
        t.eContains[i.first];
        t.eOverlaps[i.first];
        t.exonToGene[i.second.id()] = i.second.gID;
    }

    /*
     * Initialize statistics for introns
     */
    
    for (const auto &i : t.iInters.data())
    {
        t.iContains[i.first];
        t.iOverlaps[i.first];
        t.intronToGene[i.second.id()] = i.second.gID;
    }
    
    /*
     * Initialize base statistics
     */

    for (const auto &i : t.overB.hist)
    {
        t.lFPS[i.first];
        t.rFPS[i.first];
    }

    assert(!t.lFPS.empty()      && !t.rFPS.empty());
    assert(!t.eContains.empty() && !t.eOverlaps.empty());
    assert(!t.iContains.empty() && !t.iOverlaps.empty());
}

static RAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;

    RAlign::Stats stats;

    for (const auto &i : r.histGene())
    {
        initT(i.first, stats.data[i.first]);
    }

    assert(!stats.data.empty());
    return stats;
}

template <typename T> const T * matchT(const Alignment &align,
                                       Intervals<T> &inters,
                                       std::map<std::string, Counts> &contains,
                                       std::map<std::string, Counts> &overlaps,
                                       RAlign::FPStats *lFPS = nullptr,
                                       RAlign::FPStats *rFPS = nullptr)
{
    std::vector<T *> oMatches, cMatches;

    if (align.l.start == 2554775)
    {
        oMatches = oMatches;
    }
    
    /*
     * Quite likely there'll be multiple matches. Note that it's impossible to distinguish the
     * individuals due to alternative splicing. Thus, we simply increment all the possible matches.
     * Consequently, it's not possible to detect anything at the isoform level.
     */

    if (inters.contains(align.l, &cMatches))
    {
        for (const auto &i : cMatches)
        {
            contains.at(i->id())++;
        }
    }
    else
    {
        if (inters.overlap(align.l, &oMatches))
        {
            for (const auto &i : oMatches)
            {
                overlaps.at(i->id())++;
            }
        }
    }
    
    auto matches = !cMatches.empty() ? &cMatches : &oMatches;
    
    if (!matches->empty())
    {
        Base lp, rp;

        for (long i = ((*matches).size()-1); i >= 0; i--)
        {
            // Anything that fails to being mapped is counted as FP
            (*matches)[i]->map(align.l, &lp, &rp);
        }

        if (lFPS && rFPS)
        {
            lFPS->at((*matches)[0]->gID) = std::max(lFPS->at((*matches)[0]->gID), lp);
            rFPS->at((*matches)[0]->gID) = std::max(rFPS->at((*matches)[0]->gID), rp);
        }
    }
    else
    {
        /*
         * What to do for no matches? The read is outside of the reference regions and thus
         * we wouldn't have information to deduce whether this is a TP or FP. We should ignore it.
         */
    }
    
    return !cMatches.empty() ? cMatches[0] : nullptr;
}

template <typename T> void collect(const ChrID &cID,
                                   T &t,
                                   const RAlign::FPStats &lFPS,
                                   const RAlign::FPStats &rFPS,
                                   const RAlign::Options &o)
{
    /*
     * Calculating alignment statistics
     */
    
    auto aligns = [](std::map<GeneID, RAlign::MergedConfusion> &gene,
                     RAlign::MergedConfusion &over,
                     SequinHist &hist,
                     Counts n_overlaps,
                     const BinCounts &contains,
                     const BinCounts &overlaps,
                     const std::map<BinID, GeneID> &m)
    {
        /*
         * Every containment is counted as a TP
         */
        
        for (const auto &i : contains)
        {
            hist.at(m.at(i.first)) += i.second;
            gene.at(m.at(i.first)).aTP += i.second;
            over.aTP += i.second;
        }
        
        /*
         * Every overlapping is counted as a FP
         */
        
        for (const auto &i : overlaps)
        {
            gene.at(m.at(i.first)).aFP += i.second;
            over.aFP += i.second;
        }
        
        over.aFP += n_overlaps;
    };
    
    aligns(t.geneE,
           t.overE,
           t.histE,
           t.unknowns.size(),
           t.eContains,
           t.eOverlaps,
           t.exonToGene);
    
    aligns(t.geneI,
           t.overI,
           t.histI,
           0,
           t.iContains,
           t.iOverlaps,
           t.intronToGene);
    
    /*
     * Calculating statistics for each sequin (at the gene level due to alternative splicing)
     */
    
    o.info("Calculating statistics for sequins");
    
    auto genes = [](std::map<GeneID, RAlign::MergedConfusion> &gene,
                    RAlign::MergedConfusion &over,
                    const std::map<std::string, Counts> &contains,
                    const std::map<std::string, Counts> &overlaps,
                    const std::map<std::string, std::string> &m)
    {
        /*
         * Let's count number of exon/intron bins
         */
        
        for (auto &i : gene)
        {
            for (const auto &j : m)
            {
                if (i.first == j.second)
                {
                    i.second.lNR++;
                    over.lNR++;
                }
            }
        }
        
        /*
         * Every containment is counted as a TP.
         */
        
        for (const auto &i : contains)
        {
            if (i.second)
            {
                gene.at(m.at(i.first)).lTP++;
                over.lTP++;
            }
        }
    };
    
    // Do it at the exon level
    genes(t.geneE,
          t.overE,
          t.eContains,
          t.eOverlaps,
          t.exonToGene);

    // Repeat at the intron level
    genes(t.geneI,
          t.overI,
          t.iContains,
          t.iOverlaps,
          t.intronToGene);

    /*
     * Calculating statistics at the base level
     */
    
    o.info("Calculating base statistics");
    
    for (const auto &i : t.eInters.data())
    {
        auto &m  = t.geneB.at(i.second.gID);
        auto &in = i.second;
        
        const auto &gID = i.second.gID;
        
        // Update the FP at the gene level
        m.fp() = lFPS.at(gID) + rFPS.at(gID);
        
        // Update the overall FP
        t.overB.m.fp() += m.fp();
        
        Base covered = 0;
        
        in.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                covered += j - i;
                
                // Update the overall performance
                t.overB.m.tp() += j - i;
                
                // Update the distribution
                t.overB.hist.at(gID)++;
            }
        });

        m.tp() += covered;
        m.nr() += in.l().length();
        m.nq()  = m.tp() + m.fp();
        
        assert(m.nr() >= m.tp());
        
        t.overB.m.nr() += in.l().length();
        t.overB.m.nq()  = t.overB.m.tp() + t.overB.m.fp();
    }
    
    o.logInfo("Base (TP): " + std::to_string(t.overB.m.tp()));
    o.logInfo("Base (FP): " + std::to_string(t.overB.m.fp()));
    
    /*
     * Calculating for missing statistics
     */

    o.info("Calculating missing statistics");

    auto missing = [&](std::set<Missing> &misses, const BinCounts &bins)
    {
        for (const auto &bin : bins)
        {
            if (!bin.second)
            {
                misses.insert(bin.first);
            }
        }
    };

    // An intron is missing if nothing aligns to it
    missing(t.missI, t.iContains);

    // An exon is missing if nothing aligns to it
    missing(t.missE, t.eContains);
    
    /*
     * A gene is considered missing if not all it's exons have alignment
     *
     *   TODO: Need to improve the performance...
     */
    
    for (const auto &gene : t.histE)
    {
        bool missing = false;
        
        for (const auto &bin : t.eContains)
        {
            if (gene.first == t.exonToGene.at(bin.first))
            {
                if ((missing = (bin.second == 0)))
                {
                    break;
                }
            }
        }
        
        if (missing)
        {
            t.missG.insert(Missing(gene.first));
        }
    }
}

RAlign::Stats calculate(const RAlign::Options &o, Functor cal)
{
    auto stats = init();
    
    // Parsing input files
    cal(stats);

    o.info("Collecting statistics");
    
    /*
     * Collecting statistics
     */

    for (auto &i : stats.data)
    {
        collect(i.first, i.second, i.second.lFPS, i.second.rFPS, o);
    }

    /*
     * Collecting statistics for synthetic and genome
     */

    // Aggregating over all genomic chromosomes
    RAlign::MergedConfusion g_em, g_im;
    
    // Aggregating over all genomic chromosomes
    Confusion g_bm;
    
    for (auto &i : stats.data)
    {
        if (Standard::isSynthetic(i.first))
        {
            stats.s_esn = stats.sn(i.first, Metrics::AlignExon);
            stats.s_epc = stats.pc(i.first, Metrics::AlignExon);
            stats.s_isn = stats.sn(i.first, Metrics::AlignIntron);
            stats.s_ipc = stats.pc(i.first, Metrics::AlignIntron);
            
            
            
            
            std::cout << stats.data.at(i.first).eInters.size() << std::endl;
            std::cout << stats.data.at(i.first).eInters.stats().covered() << std::endl;
            
            for (const auto &x : stats.data.at(i.first).eInters.data())
            {
                x.second.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
                                                        {
                                                            if (!depth)
                                                            {
                                                                if (id == "chrT_R2_7_R2_7_2_2554764_2554890")
                                                                {
                                                                    std::cout << id << " NO " << i << "-" << j << std::endl;
                                                                    
                                                                }
                                                                
                                                                std::cout << id << " NO " << i << "-" << j << std::endl;
                                                            }
                                                            else
                                                            {
                                                                //std::cout << id << " YES " << i << "-" << j << std::endl;
                                                            }
                                                        });
            }

            
            

            
            
            stats.s_bsn = stats.sn(i.first, Metrics::AlignBase);
            stats.s_bpc = stats.pc(i.first, Metrics::AlignBase);
            stats.s_ems = stats.missProp(ChrT, MissingMetrics::MissingExon);
            stats.s_ims = stats.missProp(ChrT, MissingMetrics::MissingIntron);
            stats.s_gms = stats.missProp(ChrT, MissingMetrics::MissingGene);
            
            /*
             * Mapping from sequins to reads
             */
            
            for (const auto &j : i.second.histE)
            {
                stats.s2r[j.first] = j.second;
            }

            assert(!stats.s2r.empty());
        }
        else
        {
            const auto ems = stats.countMiss(i.first, MissingMetrics::MissingExon);
            const auto ims = stats.countMiss(i.first, MissingMetrics::MissingIntron);
            const auto gms = stats.countMiss(i.first, MissingMetrics::MissingGene);

            stats.g_ems.i += ems.i;
            stats.g_ems.n += ems.n;
            stats.g_ims.i += ims.i;
            stats.g_ims.n += ims.n;
            stats.g_gms.i += gms.i;
            stats.g_gms.n += gms.n;
            
            g_em.aTP += i.second.overE.aTP;
            g_em.aFP += i.second.overE.aFP;
            g_em.lTP += i.second.overE.lTP;
            g_em.lNR += i.second.overE.lNR;

            g_im.aTP += i.second.overI.aTP;
            g_im.aFP += i.second.overI.aFP;
            g_im.lTP += i.second.overI.lTP;
            g_im.lNR += i.second.overI.lNR;
            
            g_bm.nr() += i.second.overB.m.nr();
            g_bm.nq() += i.second.overB.m.nq();
            g_bm.tp() += i.second.overB.m.tp();
            g_bm.fp() += i.second.overB.m.fp();
            g_bm.tn() += i.second.overB.m.tn();
            g_bm.fn() += i.second.overB.m.fn();
        }
    }
    
    if (stats.data.size() > 1)
    {
        stats.g_esn = g_em.sn();
        stats.g_epc = g_em.pc();
        stats.g_isn = g_im.sn();
        stats.g_ipc = g_im.pc();
        stats.g_bsn = g_bm.sn();
        stats.g_bpc = g_bm.pc();
    }

    return stats;
}

template <typename T> const Interval * matchAlign(T &t, const Alignment &align)
{
    const Interval *match = nullptr;
    
    if (!align.spliced)
    {
        match = matchT(align,
                       t.eInters,
                       t.eContains,
                       t.eOverlaps,
                       &(t.lFPS),
                       &(t.rFPS));
    }
    else
    {
        match = matchT(align,
                       t.iInters,
                       t.iContains,
                       t.iOverlaps);
    }

    return match;
}

static void matchAlign(RAlign::Stats::Data &t,
                       const Alignment &align,
                       const ParserSAM::Info &info,
                       const RAlign::Options &o)
{
    if (!align.i && info.p.i && !(info.p.i % 1000000)) { o.wait(std::to_string(info.p.i)); }

    if (!matchAlign(t, align))
    {
        t.unknowns.push_back(UnknownAlignment(align.name, align.l));
    }
}

RAlign::Stats RAlign::analyze(const std::vector<Alignment> &aligns, const Options &o)
{
    return calculate(o, [&](RAlign::Stats &stats)
    {
        ParserSAM::Info info;
        
        for (const auto &align : aligns)
        {
            stats.update(align);

            if (!align.mapped)
            {
                return;
            }
            else if (Standard::isSynthetic(align.cID))
            {
                matchAlign(stats.data.at(align.cID), align, info, o);
            }
            else if (Standard::isGenomic(align.cID))
            {
                matchAlign(stats.data.at(align.cID), align, info, o);
            }
        }
    });
}

static void classifySyn(RAlign::Stats::Data &x,
                        const Alignment &align,
                        const ParserSAM::Info &info,
                        const RAlign::Options &o)
{
    matchAlign(x, align, info, o);
}

static void classifyGen(RAlign::Stats::Data &x,
                        const Alignment &align,
                        const ParserSAM::Info &info,
                        const RAlign::Options &o)
{
    matchAlign(x, align, info, o);
}

RAlign::Stats RAlign::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](RAlign::Stats &stats)
    {
        ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::Info &info)
        {
            stats.update(align);

            if (!align.mapped)
            {
                return;
            }            
            else if (Standard::isSynthetic(align.cID))
            {
                classifySyn(stats.data.at(align.cID), align, info, o);
            }
            else if (Standard::isGenomic(align.cID))
            {
                classifyGen(stats.data.at(align.cID), align, info, o);
            }

            if (!align.i)
            {
                if (info.spliced)
                {
                    stats.data[align.cID].n_spliced++;
                }
                else
                {
                    stats.data[align.cID].n_normal++;
                }
            }
        });
    });
}

template <typename F> std::string check(const RAlign::Stats &stats, F f, const ChrID &cID)
{
    return stats.data.count(cID) ? std::to_string(f(cID)) : "-";
}

static Scripts summary()
{
    return "-------RnaAlign Summary Statistics\n"
           "       Input alignment file: %1%\n"
           "       Reference annotation file: %2%\n\n"
           "-------Number of alignments mapped to the Synthetic and Genome\n\n"
           "       Synthetic: %3%\n"
           "       Genome:    %4%\n"
           "       Dilution:  %5$.2f\n"
           "       Unmapped:  %6%\n\n"
           "-------Reference annotation (Synthetic)\n\n"
           "       Synthetic: %7% exons\n"
           "       Synthetic: %8% introns\n"
           "       Synthetic: %9% bases\n\n"
           "-------Reference annotation (Genome)\n\n"
           "       Genome: %10% exons\n"
           "       Genome: %11% introns\n"
           "       Genome: %12% bases\n\n"
           "-------Alignments\n\n"
           "       Non-spliced (Synthetic): %13%\n"
           "       Spliced (Synthetic):     %14%\n"
           "       Total (Synthetic):       %15%\n\n"
           "       Non-spliced (Genome):    %16%\n"
           "       Spliced (Genome):        %17%\n"
           "       Total (Genome):          %18%\n\n"
           "-------Comparison of alignments to annotations (Synthetic)\n\n"
           "       *Exon level\n"
           "        Sensitivity: %19$.2f\n"
           "        Precision:   %20$.2f\n\n"
           "       *Intron level\n"
           "        Sensitivity: %21$.2f\n"
           "        Precision:   %22$.2f\n\n"
           "       *Base level\n\n"
           "        Sensitivity: %23$.2f\n"
           "        Precision:   %24$.2f\n\n"
           "       *Undetected\n"
           "        Exon:   %25$.4f\n"
           "        Intron: %26$.4f\n"
           "        Gene:   %27$.4f\n\n"
           "-------Comparison of alignments to annotations (Genome)\n\n"
           "       *Exon level\n"
           "        Sensitivity: %28%\n"
           "        Precision:   %29%\n\n"
           "       *Intron level\n"
           "        Sensitivity: %30%\n"
           "        Precision:   %31%\n\n"
           "       *Base level\n"
           "        Sensitivity: %32%\n"
           "        Precision:   %33%\n\n"
           "       *Undetected\n"
           "        Exon:   %34$.4f\n"
           "        Intron: %35$.4f\n"
           "        Gene:   %36$.4f\n";
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
    o.writer->write((boost::format(summary()) % file             // 1
                                              % GTFRef()         // 2
                                              % stats.n_syn      // 3
                                              % stats.n_gen      // 4
                                              % stats.dilution() // 5
                                              % stats.n_unmap    // 6
                                              % r.countExonSyn() // 7
                                              % r.countIntrSyn() // 8
                                              % r.countLenSyn()  // 9
                                              % CHECK(r.countExonGen()) // 10
                                              % CHECK(r.countIntrGen()) // 11
                                              % CHECK(r.countLenGen())  // 12
                                              % stats.countSpliceSyn()  // 13
                                              % stats.countNormalSyn()  // 14
                                              % (stats.countSpliceSyn() + stats.countNormalSyn())
                                              % CHECK(stats.countSpliceGen())
                                              % CHECK(stats.countNormalGen())
                                              % CHECK(stats.countSpliceGen() + stats.countNormalGen())
                                              % stats.s_esn        // 19
                                              % stats.s_epc        // 20
                                              % stats.s_isn        // 21
                                              % stats.s_ipc        // 22
                                              % stats.s_bsn        // 23
                                              % stats.s_bpc        // 24
                                              % stats.s_ems        // 25
                                              % stats.s_ims        // 26
                                              % stats.s_gms        // 27
                                              % CHECK(stats.g_esn) // 28
                                              % CHECK(stats.g_epc) // 29
                                              % CHECK(stats.g_isn) // 30
                                              % CHECK(stats.g_ipc) // 31
                                              % CHECK(stats.g_bsn) // 32
                                              % CHECK(stats.g_bpc) // 33
                                              % CHECK(stats.g_ems.percent()) // 34
                                              % CHECK(stats.g_ims.percent()) // 35
                                              % CHECK(stats.g_gms.percent()) // 36
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
    o.writer->write((boost::format(format) % "Seq"
                                           % "Reads"
                                           % "Sn_exon"
                                           % "Pc_exon"
                                           % "Sn_intron"
                                           % "Pc_intron"
                                           % "Sn_base"
                                           % "Pc_base").str());

    for (const auto &i : data.overB.hist)
    {
        // Eg: R1_1
        const auto &id = i.first;
        
        const auto &mb = data.geneB.at(id);
        const auto &me = data.geneE.at(id);
        const auto &mi = data.geneI.at(id);

        const auto m = r.findGene(ChrT, id);
        const auto reads = stats.s2r.at(id);
        
        assert(m);
        
        // Not all sequins have an intron...
        if (mi.lNR)
        {
            o.writer->write((boost::format(format) % id
                                                   % reads
                                                   % me.sn()
                                                   % me.pc()
                                                   % mi.sn()
                                                   % mi.pc()
                                                   % mb.sn()
                                                   % mb.pc()).str());
        }
        else
        {
            o.writer->write((boost::format(format) % id
                                                   % reads
                                                   % me.sn()
                                                   % me.pc()
                                                   % "--"
                                                   % "--"
                                                   % mb.sn()
                                                   % mb.pc()).str());
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