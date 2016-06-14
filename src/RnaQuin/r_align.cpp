#include "RnaQuin/r_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;
using namespace std::placeholders;

// Defined in resources.cpp
extern Scripts PlotScatter();

// Internal implementation
typedef std::function<void (RAlign::Stats &)> Functor;

typedef RAlign::Stats::AlignMetrics   AlignMetrics;
typedef RAlign::Stats::MissingMetrics MissingMetrics;

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

static RAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;

    RAlign::Stats stats;

    for (const auto &i : r.histGene())
    {
        initT(i.first, stats.data[i.first]);
    }

    assert(!stats.data.empty());
    assert(stats.data.count(ChrT));
    
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

static bool matchAlign(RAlign::Stats::Data &t,
                      const Alignment &align,
                      const ParserSAM::Info &info,
                      const RAlign::Options &o)
{
    #define REPORT_STATUS() if (!align.i && !(info.p.i % 1000000)) { o.wait(std::to_string(info.p.i)); }
    REPORT_STATUS();
    
    if (!matchAlign(t, align))
    {
        t.unknowns.push_back(UnknownAlignment(align.name, align.l));
    }
    
    return true;
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
            else
            {
                matchAlign(stats.data.at(__gID__), align, info, o);
            }
        }
    });
}

static bool classifySyn(RAlign::Stats::Data &x,
                        const Alignment &align,
                        const ParserSAM::Info &info,
                        const RAlign::Options &o)
{
    return matchAlign(x, align, info, o);
}

static bool classifyGen(RAlign::Stats::Data &x,
                        const Alignment &align,
                        const ParserSAM::Info &info,
                        const RAlign::Options &o)
{
    return Standard::isGenomic(align.cID) ? matchAlign(x, align, info, o) : false;
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
                    stats.data[align.cID].n_nspliced++;
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
           "       Dilution:  %5%\n"
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
           "        Sensitivity: %19%\n"
           "        Precision:   %20%\n\n"
           "       *Intron level\n"
           "        Sensitivity: %21%\n"
           "        Precision:   %22%\n\n"
           "       *Base level\n\n"
           "        Sensitivity: %23%\n"
           "        Precision:   %24%\n\n"
           "       *Undetected\n"
           "        Exon:   %25%\n"
           "        Intron: %26%\n"
           "        Gene:   %27%\n\n"
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
           "        Exon:   %34%\n"
           "        Intron: %35%\n"
           "        Gene:   %36%\n";
    
    
//    return "Summary for input: %1%\n\n"
//           "   ***\n"
//           "   *** Number of alignments mapped to the synthetic and genome\n"
//           "   ***\n\n"
//           "   Unmapped:  %2%\n"
//           "   Synthetic: %3% (%4%%%)\n"
//           "   Genome:    %5% (%6%%%)\n"
//           "   Dilution:  %7%\n\n"
//           "   ***\n"
//           "   *** Reference annotation (Synthetic)\n"
//           "   ***\n\n"
//           "   File: %8%\n\n"
//           "   Synthetic: %9% exons\n"
//           "   Synthetic: %10% introns\n"
//           "   Synthetic: %11% bases\n\n"
//           "   ***\n"
//           "   *** Reference annotation (Genome)\n"
//           "   ***\n\n"
//           "   File: %12%\n\n"
//           "   Genome: %13% exons\n"
//           "   Genome: %14% introns\n"
//           "   Genome: %15% bases\n\n"
//           "   ***\n"
//           "   *** Alignments\n"
//           "   ***\n\n"
//           "   Non-spliced (Synthetic): %16%\n"
//           "   Spliced (Synthetic):     %17%\n\n"
//           "   Non-spliced (Genome):    %19%\n"
//           "   Spliced (Genome):        %20%\n\n"
//           "   ***\n"
//           "   *** The following statistics are computed at the exon, intron and base level.\n"
//           "   ***\n\n"
//           "   ***                                     \n"
//           "   *** Comparison with synthetic annotation\n"
//           "   ***                                     \n\n"
//           "   -------------------- Exon level --------------------\n\n"
//           "   Sensitivity: %22%\n"
//           "   Precision:   %23%\n"
//           "   Detection Limit: %24% (attomol/ul) (%25%)\n\n"
//           "   -------------------- Intron level --------------------\n\n"
//           "   Sensitivity: %26%\n"
//           "   Precision:   %27%\n"
//           "   Detection Limit: %28% (attomol/ul) (%29%)\n\n"
//           "   -------------------- Base level --------------------\n\n"
//           "   Sensitivity: %30%\n"
//           "   Precision:   %31%\n"
//           "   Detection Limit: %32% (attomol/ul) (%33%)\n\n"
//           "   -------------------- Undetected --------------------\n\n"
//           "   Exon:   %34%\n"
//           "   Intron: %35%\n"
//           "   Gene:   %36%\n\n"
//           "   ***                                   \n"
//           "   *** Comparison with genomic annotation\n"
//           "   ***                                   \n\n"
//           "   -------------------- Exon level --------------------\n\n"
//           "   Sensitivity: %37%\n"
//           "   Precision:   %38%\n\n"
//           "   -------------------- Intron level --------------------\n\n"
//           "   Sensitivity: %39%\n"
//           "   Precision:   %40%\n\n"
//           "   -------------------- Base level --------------------\n\n"
//           "   Sensitivity: %41%\n"
//           "   Precision:   %42%\n\n"
//           "   -------------------- Undetected --------------------\n\n"
//           "   Exon:   %43%\n"
//           "   Intron: %44%\n"
//           "   Gene:   %45%\n\n";
}

static void generateSummary(const FileName &file,
                            const FileName &src,
                            const RAlign::Stats &stats,
                            const RAlign::Options &o)
{
    typedef RAlign::Stats Stats;

    const auto &r = Standard::instance().r_trans;

    #define BIND_R(x,y)   check(stats, std::bind(&x, &r, _1), y)
    #define BIND_Q(x,y)   check(stats, std::bind(&x, &stats, _1), y)
    #define BIND_E(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChrID &, enum Stats::AlignMetrics) const>(&x), &stats, _1, y), z)
    #define BIND_M(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChrID &, enum Stats::MissingMetrics) const>(&x), &stats, _1, y), z)

    const auto hasGeno = stats.data.size() > 1;
    
    o.writer->open(file);
    o.writer->write((boost::format(summary()) % "" // 1
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % "" // 10
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % "" // 20
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % ""
                                              % "" // 30
                                              % "" // 31
                                              % "" // 32
                                              % "" // 33
                                              % "" // 34
                                              % "" // 35
                                              % "").str());
    
//    o.writer->write((boost::format(summary())
//                                          % src
//                                          % stats.n_unmap
//                                          % stats.n_syn
//                                          % (100.0 * stats.synProp())
//                                          % stats.n_gen
//                                          % (100.0 * stats.genProp())
//                                          % stats.dilution()                                                // 7
//                                          % o.rAnnot                                                         // 8
//                                          % BIND_R(TransRef::countExons, ChrT)                              // 9
//                                          % BIND_R(TransRef::countIntrons, ChrT)                            // 10
//                                          % BIND_R(TransRef::exonBase, ChrT)                                // 11
//                                          % (!hasGeno ? "-"  : o.rAnnot)                                     // 12
//                                          % (!hasGeno ? "-" : BIND_R(TransRef::countExons, __gID__))        // 13
//                                          % (!hasGeno ? "-" : BIND_R(TransRef::countIntrons, __gID__))      // 14
//                                          % (!hasGeno ? "-" : BIND_R(TransRef::exonBase, __gID__))          // 15
//                                          % BIND_Q(Stats::countNSpliced, ChrT)                              // 16
//                                          % BIND_Q(Stats::countSpliced,  ChrT)                              // 17
//                                          % "" //BIND_Q(Stats::countQBases,   ChrT)                              // 18
//                                          % BIND_Q(Stats::countNSpliced, __gID__)                           // 19
//                                          % BIND_Q(Stats::countSpliced,  __gID__)                           // 20
//                                          % "" //BIND_Q(Stats::countQBases,   __gID__)                           // 21
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon, ChrT)                // 22
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon, ChrT)                // 23
//                                          % stats.limit(AlignMetrics::AlignExon).abund                      // 24
//                                          % stats.limit(AlignMetrics::AlignExon).id                         // 25
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron, ChrT)              // 26
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron, ChrT)              // 27
//                                          % stats.limit(AlignMetrics::AlignIntron).abund                    // 28
//                                          % stats.limit(AlignMetrics::AlignIntron).id                       // 29
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase, ChrT)                // 30
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase, ChrT)                // 31
//                                          % stats.limit(AlignMetrics::AlignBase).abund                      // 32
//                                          % stats.limit(AlignMetrics::AlignBase).id                         // 33
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingExon, ChrT)      // 34
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingIntron, ChrT)    // 35
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingGene, ChrT)      // 36
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon, __gID__)             // 37
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon, __gID__)             // 38
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron, __gID__)           // 39
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron, __gID__)           // 40
//                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase, __gID__)             // 41
//                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase, __gID__)             // 42
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingExon, __gID__)   // 43
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingIntron, __gID__) // 44
//                                          % BIND_M(Stats::missProp, MissingMetrics::MissingGene, __gID__)   // 45
//                     ).str());
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
    
    o.analyze("RnaAlign_quins.csv");
    writeQuins("RnaAlign_quins.csv", file, stats, o);

    /*
     * Generating RnaAlign_report.pdf
     */
    
    o.report->open("RnaAlign_report.pdf");
    o.report->addTitle("RnaAlign");
    o.report->addFile("RnaAlign_summary.stats");
    o.report->addFile("RnaAlign_quins.csv");
}