#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;
using namespace std::placeholders;

// Internal implementation
typedef std::function<void (TAlign::Stats &)> Functor;

typedef TAlign::Stats::AlignMetrics   AlignMetrics;
typedef TAlign::Stats::MissingMetrics MissingMetrics;

/*
 * -------------------- Initalization --------------------
 */

template <typename T> void initT(const ChromoID &cID, T &t)
{
    const auto &r = Standard::instance().r_trans;

    /*
     * 1. Create the structure and initalize the genes, it's different depends on the context
     */
    
    // Initalize the distributions
    t.overB.h = t.histE = t.histI = r.geneHist(cID);
    
    assert(!t.histE.empty());
    assert(!t.histI.empty());
    assert(!t.overB.h.empty());
    
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
     * Initialize statistics for exons
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

    for (const auto &i : t.overB.h)
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
    TAlign::Stats stats;

    for (const auto &cID : Standard::instance().r_trans.chromoIDs())
    {
        initT(cID, stats.data[cID]);
    }

    assert(stats.data.size() >= 1);
    
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

template <typename T> void collect(const ChromoID &cID,
                                   T &t,
                                   const TAlign::FPStats &lFPS,
                                   const TAlign::FPStats &rFPS,
                                   const TAlign::Options &o)
{
    /*
     * 1. Calculating alignment statistics.
     */
    
    o.info("Calculating alignment statistics");
    
    auto aligns = [](std::map<GeneID, TAlign::MergedConfusion> &gene,
                     TAlign::MergedConfusion &over,
                     Hist &h,
                     Counts unknowns,
                     const BinCounts &contains,
                     const BinCounts &overlaps,
                     const std::map<BinID, GeneID> &m)
    {
        /*
         * Every containment is counted as a TP.
         */
        
        for (const auto &i : contains)
        {
            h.at(m.at(i.first)) += i.second;
            gene.at(m.at(i.first)).aTP += i.second;
            over.aTP += i.second;
        }
        
        /*
         * Every overlapping is counted as a FP.
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
        
        // Update the FP at the overall level
        t.overB.m.fp() += m.fp();
        
        Base covered = 0;
        
        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                covered += j - i;
                
                // Update the overall performance
                t.overB.m.tp() += j - i;
                
                // Update the distribution
                t.overB.h.at(gID)++;
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
    /*
     * 1: Initalize the statistics
     */
    
    TAlign::Stats stats = init();
    
    /*
     * 2: Parsing the inputs. For instance, parsing an input file.
     */
    
    cal(stats);

    o.info("Collecting statistics");
    
    /*
     * 3: Collecting statistics
     */

    for (auto &i : stats.data)
    {
        collect(i.first, i.second, i.second.lFPS, i.second.rFPS, o);
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

static void classifyEndo(TAlign::Stats::Data &t,
                         const Alignment &align,
                         const ParserSAM::AlignmentInfo &info,
                         const TAlign::Options &o)
{
    REPORT_STATUS();
    
    if (!align.mapped)
    {
        return;
    }

    if (!matchAlign(t, align))
    {
        t.unknowns.push_back(UnknownAlignment(align.qName, align.l));
    }
}

static void classifyChrT(TAlign::Stats::Data &t,
                         const Alignment &align,
                         const ParserSAM::AlignmentInfo &info,
                         const TAlign::Options &o)
{
    assert(align.id == Standard::chrT);
    
    REPORT_STATUS();
    
    if (!align.mapped)
    {
        return;
    }
    else if (!matchAlign(t, align))
    {
        t.unknowns.push_back(UnknownAlignment(align.qName, align.l));
    }
}

TAlign::Stats TAlign::analyze(const std::vector<Alignment> &aligns, const Options &o)
{
    return calculate(o, [&](TAlign::Stats &stats)
    {
        ParserSAM::AlignmentInfo info;
        
        for (const auto &align : aligns)
        {
            stats.update(align);

            if (align.id == ChrT)
            {
                classifyChrT(stats.data.at(ChrT), align, info, o);
            }
            else
            {
                classifyEndo(stats.data.at(align.id), align, info, o);
            }
        }
    });
}

TAlign::Stats TAlign::analyze(const FileName &file, const Options &o)
{
    // We'll need the factors for generation.
    assert(o.exp);
    
    o.analyze(file);
    
    return calculate(o, [&](TAlign::Stats &stats)
    {
        ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
        {
            stats.update(align);

            if (align.id == ChrT)
            {
                classifyChrT(stats.data.at(ChrT), align, info, o);
            }
            else if (stats.data.count(align.id))
            {
                classifyEndo(stats.data.at(align.id), align, info, o);
            }
            
            /*
             * Any read that is not aligned into the reference annoation is worthless. We don't know if the locus is an exon
             * or an intron... We can't really do much...
             */
        });
    });
}

template <typename F> std::string check(const TAlign::Stats &stats, F f, const ChromoID &cID)
{
    const auto id = cID == ChrT ? ChrT : "chr1"; // TODO: Fix this...
    
    // Don't try unless we're sure it'll work
    if (!stats.data.count(id))
    {
        return "????";
    }

    return std::to_string(f(id));
}

/*
 * Summary files for the replicates and samples
 */

static Scripts replicateSummary()
{
    return "Summary for input: %1%\n\n"
           "   ***\n"
           "   *** Proportion of reads mapped to the synthetic and experiment\n"
           "   ***\n\n"
           "   Unmapped:   %2% reads\n"
           "   Synthetic:  %3% (%4%%%) reads\n"
           "   Experiment: %5% (%6%%%) reads\n"
           "   Dilution:   %7%\n\n"
           "   ***\n"
           "   *** Reference annotation (Synthetic)\n"
           "   ***\n\n"
           "   File: %8%\n\n"
           "   Synthetic: %9% exons\n"
           "   Synthetic: %10% introns\n"
           "   Synthetic: %11% bases\n\n"
           "   ***\n"
           "   *** Reference annotation (Experiment)\n"
           "   ***\n\n"
           "   File: %12%\n\n"
           "   Experiment: %13% exons\n"
           "   Experiment: %14% introns\n"
           "   Experiment: %15% bases\n\n"
           "   ***\n"
           "   *** Input file\n"
           "   ***\n\n"
           "   Synthetic:  %16% exons\n"
           "   Synthetic:  %17% introns\n"
           "   Synthetic:  %18% bases\n\n"
           "   Experiment: %19% exons\n"
           "   Experiment: %20% introns\n"
           "   Experiment: %21% bases\n\n"
           "   ***\n"
           "   *** The following statistics are computed at the exon, intron and base level.\n"
           "   ***\n"
           "   *** Exon level is defined by performance per exon. An alignment that\n"
           "   *** is not mapped entirely within an exon is considered as a FP. The\n"
           "   *** intron level is similar.\n"
           "   ***\n"
           "   *** Base level is defined by performance per nucleotide. A partial\n"
           "   *** mapped read will have FP and TP.\n"
           "   ***\n"
           "   *** Please refer to the paper \"Evaluation of gene structure prediction programs\" for more details\n"
           "   ***\n\n"
           "   *************************************************\n"
           "   ***                                           ***\n"
           "   ***    Statistics for synthetic chromosome    ***\n"
           "   ***                                           ***\n"
           "   *************************************************\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %22%\n"
           "   Specificity: %23%\n"
           "   Detection Limit: %24% (attomol/ul) (%25%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %26%\n"
           "   Specificity: %27%\n"
           "   Detection Limit: %28% (attomol/ul) (%29%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %30%\n"
           "   Specificity: %31%\n"
           "   Detection Limit: %32% (attomol/ul) (%33%)\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %34%\n"
           "   Intron: %35%\n"
           "   Gene:   %36%\n\n"
           "   ***************************************\n"
           "   ***                                 ***\n"
           "   ***    Statistics for experiment    ***\n"
           "   ***                                 ***\n"
           "   ***************************************\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %37%\n"
           "   Specificity: %38%\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %39%\n"
           "   Specificity: %40%\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %41%\n"
           "   Specificity: %42%\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %43%\n"
           "   Intron: %44%\n"
           "   Gene:   %45%\n\n";
}

static void writeSummary(const FileName &file, const FileName &src, const TAlign::Stats &stats, const TAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    typedef TAlign::Stats Stats;

    #define BIND_R(x,y)   check(stats, std::bind(&x, &r, _1), y)
    #define BIND_Q(x,y)   check(stats, std::bind(&x, &stats, _1), y)
    #define BIND_E(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChromoID &, enum Stats::AlignMetrics) const>(&x), &stats, _1, y), z)
    #define BIND_M(x,y,z) check(stats, std::bind(static_cast<double (Stats::*)(const ChromoID &, enum Stats::MissingMetrics) const>(&x), &stats, _1, y), z)

    o.writer->open(file);
    o.writer->write((boost::format(replicateSummary())
                                          % src
                                          % stats.unmapped
                                          % stats.n_chrT
                                          % (100.0 * stats.chrTProp())
                                          % stats.n_endo
                                          % (100.0 * stats.endoProp())
                                          % stats.dilution()                                                // 7
                                          % extractFile(o.rChrT)                                            // 8
                                          % BIND_R(TransRef::countExons, ChrT)                              // 9
                                          % BIND_R(TransRef::countIntrons, ChrT)                            // 10
                                          % BIND_R(TransRef::exonBase, ChrT)                                // 11
                                          % (o.rEndo.empty() ? "-" : extractFile(o.rEndo))                  // 12
                                          % (o.rEndo.empty() ? "NA" : BIND_R(TransRef::countExons, Endo))   // 13
                                          % (o.rEndo.empty() ? "NA" : BIND_R(TransRef::countIntrons, Endo)) // 14
                                          % (o.rEndo.empty() ? "NA" : BIND_R(TransRef::exonBase, Endo))     // 15
                                          % BIND_Q(Stats::qExons, ChrT)                                     // 16
                                          % BIND_Q(Stats::qIntrons, ChrT)                                   // 17
                                          % BIND_Q(Stats::qBases, ChrT)                                     // 18
                                          % BIND_Q(Stats::qExons, Endo)                                     // 19
                                          % BIND_Q(Stats::qIntrons, Endo)                                   // 20
                                          % BIND_Q(Stats::qBases, Endo)                                     // 21
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
                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon, Endo)                // 37
                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon, Endo)                // 38
                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron, Endo)              // 39
                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron, Endo)              // 40
                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase, Endo)                // 41
                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase, Endo)                // 42
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingExon, Endo)      // 43
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingIntron, Endo)    // 44
                                          % BIND_M(Stats::missProp, MissingMetrics::MissingGene, Endo)      // 45
                     ).str());
    o.writer->close();
}

static void writeSequins(const FileName &file, const FileName &src, const TAlign::Stats &stats, const TAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    o.writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->write((boost::format(format) % "ID"
                                           % "Abundance (attomol/ul)"
                                           % "Covered"
                                           % "Sensitivity (Exon)"
                                           % "Specificity (Exon)"
                                           % "Sensitivity (Intron)"
                                           % "Specificity (Intron)"
                                           % "Sensitivity (Base)"
                                           % "Specificity (Base)").str());

    for (const auto &i : stats.data.at(ChrT).overB.h)
    {
        Base length   = 0;
        Base nonZeros = 0;
        
        for (const auto &j : stats.data.at(ChrT).eInters.data())
        {
            if (j.second.gID == i.first)
            {
                const auto eStats = j.second.stats();
                
                length   += eStats.length;
                nonZeros += eStats.nonZeros;
                
                assert(length >= nonZeros);
            }
        }
        
        const auto covered = static_cast<double>(nonZeros) / length;
        
        const auto &mb = stats.data.at(ChrT).geneB.at(i.first);
        const auto &me = stats.data.at(ChrT).geneE.at(i.first);
        const auto &mi = stats.data.at(ChrT).geneI.at(i.first);
        
        // Not all sequins have an intron...
        if (mi.lNR)
        {
            o.writer->write((boost::format(format) % i.first
                                                   % r.findGene(ChrT, i.first)->abund(Mix_1)
                                                   % covered
                                                   % me.sn()
                                                   % me.pc()
                                                   % mi.sn()
                                                   % mi.pc()
                                                   % mb.sn()
                                                   % mb.pc()).str());
        }
        else
        {
            o.writer->write((boost::format(format) % i.first
                                                   % r.findGene(ChrT, i.first)->abund(Mix_1)
                                                   % covered
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

void TAlign::report(const std::vector<FileName> &files, const Options &o)
{
    const auto stats = TAlign::analyze(files, o);
    
    for (auto i = 0; i < files.size(); i++)
    {
        const auto &name = o.exp->names().at(i);
        
        // Create the directory if haven't
        o.writer->create(name);
        
        // Generating summary statistics for the replicate
        writeSummary(name + "/TransAlign_summary.stats", extractFile(files[i]), stats[i], o);
        
        // Generating sequin statistics for the replicate
        writeSequins(name + "/TransAlign_quins.stats", extractFile(files[i]), stats[i], o);
    }
}