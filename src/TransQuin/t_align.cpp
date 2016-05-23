#include "TransQuin/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;
using namespace std::placeholders;

// Defined in resources.cpp
extern Scripts PlotTReads();

// Internal implementation
typedef std::function<void (TAlign::Stats &)> Functor;

typedef TAlign::Stats::AlignMetrics   AlignMetrics;
typedef TAlign::Stats::MissingMetrics MissingMetrics;

// Defined for convenience
static ChrID __gID__;

template <typename T> void initT(const ChrID &cID, T &t)
{
    const auto &r = Standard::instance().r_trans;

    /*
     * Create the structure and initalize the genes, it's different depends on the context
     */
    
    // Initalize the distributions
    t.overB.hist = t.histE = t.histI = r.geneHist(cID);
    
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

static TAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;

    TAlign::Stats stats;

    initT(ChrT, stats.data[ChrT]);

    if (!r.genoID().empty())
    {
        initT(r.genoID(), stats.data[__gID__ = r.genoID()]);
    }
    
    assert(!stats.data.empty());
    
    return stats;
}

template <typename T> const T * matchT(const Alignment &align,
                                       Intervals<T> &inters,
                                       std::map<std::string, Counts> &contains,
                                       std::map<std::string, Counts> &overlaps,
                                       TAlign::FPStats *lFPS = nullptr,
                                       TAlign::FPStats *rFPS = nullptr)
{
    std::vector<T *> oMatches, cMatches;

    /*
     * It's quite likely there'll be more than a match. Note that it's impossible to distinguish the
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
        
        // Anything that fails to being mapped is counted as FP
        (*matches)[0]->map(align.l, &lp, &rp);

        if (lFPS && rFPS)
        {
            lFPS->at((*matches)[0]->gID) = std::max(lFPS->at((*matches)[0]->gID), lp);
            rFPS->at((*matches)[0]->gID) = std::max(rFPS->at((*matches)[0]->gID), rp);
        }
    }
    
    return !cMatches.empty() ? cMatches[0] : nullptr;
}

template <typename T> void collect(const ChrID &cID,
                                   T &t,
                                   const TAlign::FPStats &lFPS,
                                   const TAlign::FPStats &rFPS,
                                   const TAlign::Options &o)
{
    /*
     * Calculating alignment statistics
     */
    
    auto aligns = [](std::map<GeneID, TAlign::MergedConfusion> &gene,
                     TAlign::MergedConfusion &over,
                     SequinHist &hist,
                     Counts unknowns,
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
        
        over.aFP += unknowns;
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
     * 2. Calculating statistics for each sequin (at the gene level due to alternative splicing)
     */
    
    o.info("Calculating statistics for sequins");
    
    auto genes = [](std::map<GeneID, TAlign::MergedConfusion> &gene,
                    TAlign::MergedConfusion &over,
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
     * 3. Calculating metrics at the base level.
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
    
    // An exon is missing if no alignment aligns to it
    missing(t.missE, t.eContains);
    
    // An intron is missing if no alignment aligns to it
    missing(t.missI, t.iContains);
    
    /*
     * A gene is considered missing if not all it's exons have alignment
     *
     *   TODO: Need to improve the performance...
     */
    
//    if (cID == ChrT)
    {
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
}

TAlign::Stats calculate(const TAlign::Options &o, Functor cal)
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
     * Mapping from sequins to reads
     */

    for (const auto &i : stats.data.at(ChrT).histE)
    {
        stats.s2r[i.first] = i.second;
    }
    
    assert(!stats.s2r.empty());
    
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

static bool matchAlign(TAlign::Stats::Data &t,
                      const Alignment &align,
                      const ParserSAM::AlignmentInfo &info,
                      const TAlign::Options &o)
{
    #define REPORT_STATUS() if (!align.i && !(info.p.i % 1000000)) { o.wait(std::to_string(info.p.i)); }
    REPORT_STATUS();
    
    if (!matchAlign(t, align))
    {
        t.unknowns.push_back(UnknownAlignment(align.name, align.l));
        return true;
    }
    
    return false;
}

TAlign::Stats TAlign::analyze(const std::vector<Alignment> &aligns, const Options &o)
{
    return calculate(o, [&](TAlign::Stats &stats)
    {
        ParserSAM::AlignmentInfo info;
        
        for (const auto &align : aligns)
        {
            stats.update(align);

            if (!align.mapped)
            {
                return;
            }
            else if (align.cID == ChrT)
            {
                matchAlign(stats.data.at(ChrT), align, info, o);
            }
            else
            {
                matchAlign(stats.data.at(__gID__), align, info, o);
            }
        }
    });
}

TAlign::Stats TAlign::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](TAlign::Stats &stats)
    {
        ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
        {
            stats.update(align);

            if (!align.mapped)
            {
                return;
            }
            
            ChrID cID;
            bool succeed = false;

            if (align.cID == ChrT)
            {
                succeed = matchAlign(stats.data.at(cID = ChrT), align, info, o);
            }
            else if (stats.data.count(align.cID))
            {
                succeed = matchAlign(stats.data.at(cID = align.cID), align, info, o);
            }
            
            if (succeed && !align.i)
            {
                if (info.spliced)
                {
                    stats.data[cID].spliced++;
                }
                else
                {
                    stats.data[cID].nspliced++;
                }
            }

            /*
             * Any read that is not covered by the reference annoation is worthless. We don't know if the
             * aligned position is an exon or an intron... We can't really do much...
             */
        });
    });
}

template <typename F> std::string check(const TAlign::Stats &stats, F f, const ChrID &cID)
{
    return stats.data.count(cID) ? std::to_string(f(cID)) : "-";
}

static Scripts replicateSummary()
{
    return "Summary for input: %1%\n\n"
           "   ***\n"
           "   *** Number of alignments mapped to the synthetic and genome\n"
           "   ***\n\n"
           "   Unmapped:  %2%\n"
           "   Synthetic: %3% (%4%%%)\n"
           "   Genome:    %5% (%6%%%)\n"
           "   Dilution:  %7%\n\n"
           "   ***\n"
           "   *** Reference annotation (Synthetic)\n"
           "   ***\n\n"
           "   File: %8%\n\n"
           "   Synthetic: %9%  exons\n"
           "   Synthetic: %10% introns\n"
           "   Synthetic: %11% bases\n\n"
           "   ***\n"
           "   *** Reference annotation (Genome)\n"
           "   ***\n\n"
           "   File: %12%\n\n"
           "   Genome: %13% exons\n"
           "   Genome: %14% introns\n"
           "   Genome: %15% bases\n\n"
           "   ***\n"
           "   *** Alignments\n"
           "   ***\n\n"
           "   Non-spliced (Synthetic): %16%\n"
           "   Spliced (Synthetic):     %17%\n\n"
           "   Non-spliced (Genome):    %19%\n"
           "   Spliced (Genome):        %20%\n\n"
           "   ***\n"
           "   *** The following statistics are computed at the exon, intron and base level.\n"
           "   ***\n\n"
           "   ***                                     \n"
           "   *** Comparison with synthetic annotation\n"
           "   ***                                     \n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %22%\n"
           "   Precision:   %23%\n"
           "   Detection Limit: %24% (attomol/ul) (%25%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %26%\n"
           "   Precision:   %27%\n"
           "   Detection Limit: %28% (attomol/ul) (%29%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %30%\n"
           "   Precision:   %31%\n"
           "   Detection Limit: %32% (attomol/ul) (%33%)\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %34%\n"
           "   Intron: %35%\n"
           "   Gene:   %36%\n\n"
           "   ***                                   \n"
           "   *** Comparison with genomic annotation\n"
           "   ***                                   \n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %37%\n"
           "   Precision:   %38%\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %39%\n"
           "   Precision:   %40%\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %41%\n"
           "   Precision:   %42%\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %43%\n"
           "   Intron: %44%\n"
           "   Gene:   %45%\n\n";
}

static void generateSummary(const FileName &file,
                            const FileName &src,
                            const TAlign::Stats &stats,
                            const TAlign::Options &o)
{
    typedef TAlign::Stats Stats;

    const auto &r = Standard::instance().r_trans;

    #define BIND_R(x,y)   check(stats, std::bind(&x, &r, _1), y)
    #define BIND_Q(x,y)   check(stats, std::bind(&x, &stats, _1), y)
    #define BIND_E(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChrID &, enum Stats::AlignMetrics) const>(&x), &stats, _1, y), z)
    #define BIND_M(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChrID &, enum Stats::MissingMetrics) const>(&x), &stats, _1, y), z)

    const auto hasGeno = !o.rGeno.empty();
    
    o.writer->open(file);
    o.writer->write((boost::format(replicateSummary())
                                          % src
                                          % stats.unmapped
                                          % stats.n_chrT
                                          % (100.0 * stats.chrTProp())
                                          % stats.n_geno
                                          % (100.0 * stats.endoProp())
                                          % stats.dilution()                                                // 7
                                          % o.rChrT                                                         // 8
                                          % BIND_R(TransRef::countExons, ChrT)                              // 9
                                          % BIND_R(TransRef::countIntrons, ChrT)                            // 10
                                          % BIND_R(TransRef::exonBase, ChrT)                                // 11
                                          % (!hasGeno ? "-"  : o.rGeno)                                     // 12
                                          % (!hasGeno ? "-" : BIND_R(TransRef::countExons, __gID__))        // 13
                                          % (!hasGeno ? "-" : BIND_R(TransRef::countIntrons, __gID__))      // 14
                                          % (!hasGeno ? "-" : BIND_R(TransRef::exonBase, __gID__))          // 15
                                          % BIND_Q(Stats::countNSpliced, ChrT)                              // 16
                                          % BIND_Q(Stats::countSpliced,  ChrT)                              // 17
                                          % "" //BIND_Q(Stats::countQBases,   ChrT)                              // 18
                                          % BIND_Q(Stats::countNSpliced, __gID__)                           // 19
                                          % BIND_Q(Stats::countSpliced,  __gID__)                           // 20
                                          % "" //BIND_Q(Stats::countQBases,   __gID__)                           // 21
                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon, ChrT)                // 22
                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon, ChrT)                // 23
                                          % stats.limit(AlignMetrics::AlignExon).abund                      // 24
                                          % stats.limit(AlignMetrics::AlignExon).id                         // 25
                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron, ChrT)              // 26
                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron, ChrT)              // 27
                                          % stats.limit(AlignMetrics::AlignIntron).abund                    // 28
                                          % stats.limit(AlignMetrics::AlignIntron).id                       // 29
                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase, ChrT)                // 30
                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase, ChrT)                // 31
                                          % stats.limit(AlignMetrics::AlignBase).abund                      // 32
                                          % stats.limit(AlignMetrics::AlignBase).id                         // 33
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingExon, ChrT)      // 34
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingIntron, ChrT)    // 35
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingGene, ChrT)      // 36
                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon, __gID__)             // 37
                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon, __gID__)             // 38
                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron, __gID__)           // 39
                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron, __gID__)           // 40
                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase, __gID__)             // 41
                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase, __gID__)             // 42
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingExon, __gID__)   // 43
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingIntron, __gID__) // 44
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingGene, __gID__)   // 45
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file,
                       const FileName &src,
                       const TAlign::Stats &stats,
                       const TAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";

    const auto &data = stats.data.at(ChrT);
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "seq"
                                           % "input"
                                           % "reads"
                                           % "sn_exon"
                                           % "pc_exon"
                                           % "sn_intron"
                                           % "pc_intron"
                                           % "sn_base"
                                           % "pc_base").str());

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
                                                   % m->concent(Mix_1)
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
                                                   % m->concent(Mix_1)
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

void TAlign::report(const FileName &file, const Options &o)
{
    const auto stats = TAlign::analyze(file, o);
    
    o.info("Generating statistics");
    
    /*
     * Generating TransAlign_summary.stats
     */
    
    o.analyze("TransAlign_summary.stats");
    generateSummary("TransAlign_summary.stats", file, stats, o);

    /*
     * Generating TransAlign_quins.stats
     */
    
    o.analyze("TransAlign_quins.stats");
    writeQuins("TransAlign_quins.stats", file, stats, o);

    /*
     * Generating TransAlign_reads.R
     */
    
    o.generate("TransAlign_reads.R");
    o.writer->open("TransAlign_reads.R");
    o.writer->write(RWriter::createScript("TransAlign_quins.stats", PlotTReads()));
    o.writer->close();

    /*
     * Generating TransAlign_report.pdf
     */
    
    o.report->open("TransAlign_report.pdf");
    o.report->addTitle("TransAlign");
    o.report->addFile("TransAlign_summary.stats");
    o.report->addFile("TransAlign_quins.stats");
    o.report->addFile("TransAlign_reads.R");
}