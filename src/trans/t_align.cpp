#include "trans/t_align.hpp"
#include "data/accumulator.hpp"
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

// Template function used by init()
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
     * Initialize exon statistics
     */
    
    for (const auto &i : t.eInters.data())
    {
        t.eContains[i.first];
        t.eOverlaps[i.first];
        t.exonToGene[i.second.id()] = i.second.gID;
    }

    /*
     * Initialize intron statistics
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

    /*
     * Initalize for each chromosome. The results will be pooled together.
     */
    
    for (const auto &cID : Standard::instance().r_trans.chromoIDs())
    {
        initT(cID, stats.data[cID]);
    }

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
     * individuals due to alternative splicing. Thus, we simply increment for all the possible matches.
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
    
    o.info("Base (TP): " + std::to_string(t.overB.m.tp()));
    o.info("Base (FP): " + std::to_string(t.overB.m.fp()));
    
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
     * A gene is considered missing if not all exons have alignment aligned to it... 
     *
     *   TODO: Need to improve the performance...
     */
    
    if (cID == ChrT)
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

/*
 * This function provides flexibiltiy in combining values across conditions. It can also be used for synthetic chromosome.
 */

template <typename Data, typename F> std::string combine(const Data &data, F f)
{
    // What's the value for the synthetic chromosome?
    const auto chrT = std::to_string(f(ChrT));
    
    return chrT;
}

/*
 * Summary statistics for a single replicate.
 */

static std::string replicateSummary()
{
    return "Summary for file: %1%\n\n"
           "   Unmapped:   %2% reads\n"
           "   Experiment: %3% (%4%%%) reads\n"
           "   Synthetic:  %5% (%6%%%) reads\n\n"
           "   Reference:  %7% exons\n"
           "   Reference:  %8% introns\n"
           "   Reference:  %9% bases\n\n"
           "   Query:      %10% exons\n"
           "   Query:      %11% introns\n"
           "   Query:      %12% bases\n\n"
           "   Dilution:   %13%\n\n"
           "   ***\n"
           "   *** The following statistics are computed at the exon, intron and base level.\n"
           "   ***\n"
           "   *** Exon level is defined by performance per exon. An alignment that\n"
           "   *** is not mapped entirely within an exon is considered as a FP. The\n"
           "   *** intron level is similar.\n"
           "   ***\n"
           "   *** Base level is defined by performance per nucleotide. A partial\n"
           "   *** mapped read will have FP and TP.\n"
           "   ***\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %14%\n"
           "   Specificity: %15%\n"
           "   Detection:   %16% (%17%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %18%\n"
           "   Specificity: %19%\n"
           "   Detection:   %20% (%21%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %22%\n"
           "   Specificity: %23%\n"
           "   Detection:   %24% (%25%)\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %26%\n"
           "   Intron: %27%\n"
           "   Gene:   %28%\n\n";
}

static void writeSummary(const FileName &file, const FileName &src, const TAlign::Stats &stats, std::shared_ptr<Writer> writer)
{
    const auto &r = Standard::instance().r_trans;

    typedef TAlign::Stats Stats;

    #define BIND_R(x)    combine(stats.data, std::bind(&x, &r, _1))
    #define BIND_Q(x)    combine(stats.data, std::bind(&x, &stats, _1))
    #define BIND_E(x, y) combine(stats.data, std::bind(static_cast<double (Stats::*)(const ChromoID &, enum Stats::AlignMetrics) const>(&x), &stats, _1, y))
    #define BIND_M(x, y) combine(stats.data, std::bind(static_cast<double (Stats::*)(const ChromoID &, enum Stats::MissingMetrics) const>(&x), &stats, _1, y))

    writer->open(file);
    writer->write((boost::format(replicateSummary())
                                          % src
                                          % stats.unmapped
                                          % stats.n_endo
                                          % (100.0 * stats.endoProp())
                                          % stats.n_chrT
                                          % (100.0 * stats.chrTProp())                            // 6
                                          % BIND_R(TransRef::countExons)                          // 7
                                          % BIND_R(TransRef::countIntrons)                        // 8
                                          % BIND_R(TransRef::exonBase)                            // 9
                                          % BIND_Q(Stats::qExons)                                 // 10
                                          % BIND_Q(Stats::qIntrons)                               // 11
                                          % BIND_Q(Stats::qBases)                                 // 12
                                          % stats.dilution()                                      // 13
                                          % BIND_E(Stats::sn, AlignMetrics::AlignExon)            // 14
                                          % BIND_E(Stats::pc, AlignMetrics::AlignExon)            // 15
                                          % stats.limit(AlignMetrics::AlignExon).abund            // 16
                                          % stats.limit(AlignMetrics::AlignExon).id               // 17
                                          % BIND_E(Stats::sn, AlignMetrics::AlignIntron)          // 18
                                          % BIND_E(Stats::pc, AlignMetrics::AlignIntron)          // 19
                                          % stats.limit(AlignMetrics::AlignIntron).abund          // 20
                                          % stats.limit(AlignMetrics::AlignIntron).id             // 21
                                          % BIND_E(Stats::sn, AlignMetrics::AlignBase)            // 22
                                          % BIND_E(Stats::pc, AlignMetrics::AlignBase)            // 23
                                          % stats.limit(AlignMetrics::AlignBase).abund            // 24
                                          % stats.limit(AlignMetrics::AlignBase).id               // 25
                                          % BIND_M(Stats::missPercent, MissingMetrics::MissingExon)
                                          % BIND_M(Stats::missPercent, MissingMetrics::MissingIntron)
                                          % BIND_M(Stats::missPercent, MissingMetrics::MissingGene)
                     ).str());
    writer->close();
}

static void writeSequins(const FileName &file, const FileName &src, const TAlign::Stats &stats, std::shared_ptr<Writer> writer)
{
    writer->open(file);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    
    writer->write((boost::format(format) % "ID"
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
            writer->write((boost::format(format) % i.first
                                                 % covered
                                                 % me.sn()
                                                 % me.precise()
                                                 % mi.sn()
                                                 % mi.precise()
                                                 % mb.sn()
                                                 % mb.ac()).str());
        }
        else
        {
            writer->write((boost::format(format) % i.first
                                                 % me.sn()
                                                 % me.precise()
                                                 % "--"
                                                 % "--"
                                                 % mb.sn()
                                                 % mb.ac()
                                                 % covered).str());
        }
    }
    
    writer->close();
}

/*
 * Write summary statistics for a replicate. file is the file name of the replicate. name is the name of the replicate.
 *
 *     Eg: writeReplicate(..., "A1/accepted_hits.bam", "A1", ...)
 */

static void writeReplicate(const TAlign::Stats &stats, const FileName &file, const std::string &name, const TAlign::Options &o)
{
    // Create the directory if haven't
    o.writer->create(name);
    
    // Generating summary statistics for the replicate
    writeSummary(name + "/TransAlign_summary.stats", extractFile(file), stats, o.writer);
    
    // Generating sequin statistics for the replicate
    writeSequins(name + "/TransAlign_quins.stats", extractFile(file), stats, o.writer);
}

static std::string pooledSummary()
{
    return "Summary for file: %1%\n\n"
           "   Unmapped:   %2% reads\n"
           "   Experiment: %3% (%4%%%) reads\n"
           "   Synthetic:  %5% (%6%%%) reads\n\n"
           "   Query:      %7% exons\n"
           "   Query:      %8% introns\n"
           "   Query:      %9% bases\n\n"
           "   ***\n"
           "   *** The following statistics are computed at the exon, intron and base level.\n"
           "   ***\n"
           "   *** Exon level is defined by performance per exon. An alignment that\n"
           "   *** is not mapped entirely within an exon is considered as a FP. The\n"
           "   *** intron level is similar.\n"
           "   ***\n"
           "   *** Base level is defined by performance per nucleotide. A partial\n"
           "   *** mapped read will have FP and TP.\n"
           "   ***\n\n"
           "   -------------------- Exon level --------------------\n\n"
           "   Sensitivity: %10%\n"
           "   Specificity: %11%\n"
           "   Detection:   %12% (%13%)\n\n"
           "   -------------------- Intron level --------------------\n\n"
           "   Sensitivity: %14%\n"
           "   Specificity: %15%\n"
           "   Detection:   %16% (%17%)\n\n"
           "   -------------------- Base level --------------------\n\n"
           "   Sensitivity: %18%\n"
           "   Specificity: %19%\n\n"
           "   Detection:   %20% (%21%)\n\n"
           "   -------------------- Undetected --------------------\n\n"
           "   Exon:   %22%\n"
           "   Intron: %23%\n"
           "   Gene:   %24%\n\n";
}

void TAlign::report(const std::vector<FileName> &files, const Options &o)
{
    std::map<ChromoID, Accumulator> accs;
    
    const auto stats = TAlign::analyze(files, o);
    
    /*
     * Process each replicate in orders. Later, we'll pool the information to generate a summary for all replicates.
     * In order to calculate the variation between replicates, we'll add them to an accumulator.
     */
    
    for (auto i = 0; i < files.size(); i++)
    {
        const auto &stat = stats[i];
        
        accs[ExpT].add("n_endo",    stat.n_endo);
        accs[ExpT].add("n_chrT",    stat.n_chrT);
        accs[ExpT].add("dilution",  stat.dilution());
        accs[ExpT].add("unmapped",  stat.unmapped);
        accs[ExpT].add("expTProp",  stat.endoProp());
        accs[ExpT].add("chrTProp",  stat.chrTProp());
        accs[ExpT].add("limitE",    stat.limit(AlignMetrics::AlignExon));
        accs[ExpT].add("limitI",    stat.limit(AlignMetrics::AlignIntron));

        auto f = [&](const ChromoID &id)
        {
            accs[id].add("qExons",   stat.qExons(id));
            accs[id].add("qIntrons", stat.qIntrons(id));
            accs[id].add("qBases",   stat.qBases(id));
            accs[id].add("snE",      stat.sn(id, AlignMetrics::AlignExon));
            accs[id].add("pcE",      stat.pc(id, AlignMetrics::AlignExon));
            accs[id].add("snI",      stat.sn(id, AlignMetrics::AlignIntron));
            accs[id].add("pcI",      stat.pc(id, AlignMetrics::AlignIntron));
            accs[id].add("snB",      stat.sn(id, AlignMetrics::AlignBase));
            accs[id].add("pcB",      stat.pc(id, AlignMetrics::AlignBase));
            accs[id].add("limitB",   stat.data.at(id).overB.limit);
            accs[id].add("missE",    stat.missPercent(id, MissingMetrics::MissingExon));
            accs[id].add("missI",    stat.missPercent(id, MissingMetrics::MissingIntron));
            accs[id].add("missG",    stat.missPercent(id, MissingMetrics::MissingGene));
        };

        for (const auto &cID : Standard::instance().r_trans.chromoIDs())
        {
            f(cID);
        }
        
        // Generate summary statistic for the replicate
        writeReplicate(stats[i], files[i], o.exp->names().at(i), o);
    }

    /*
     * Generating pooled summary statistics
     */

    std::string concated;
    
    for (const auto &file : files)
    {
        if (concated.empty())
        {
            concated = file;
        }
        else
        {
            concated = concated + "\n                     " + file;
        }
    }
    
    auto f = [&](const ChromoID &id)
    {
        const auto acc = accs.at(id);
        
        o.writer->write((boost::format(pooledSummary()) % concated
                                                        % acc.value("unmapped")()
                                                        % acc.value("n_endo")()
                                                        % acc.value("expTProp")()    // 4
                                                        % acc.value("n_chrT")()      // 5
                                                        % acc.value("chrTProp")()    // 6
                                                        % acc.value("qExons")()      // 7
                                                        % acc.value("qIntrons")()    // 8
                                                        % acc.value("qBases")()      // 9
                                                        % acc.value("snE")()         // 10
                                                        % acc.value("pcE")()         // 11
                                                        % acc.limits("limitE").abund // 12
                                                        % acc.limits("limitE").id    // 13
                                                        % acc.value("snI")()         // 14
                                                        % acc.value("pcI")()         // 15
                                                        % acc.limits("limitI").abund // 16
                                                        % acc.limits("limitI").id    // 17
                                                        % acc.value("snB")()         // 18
                                                        % acc.value("pcB")()         // 19
                                                        % acc.limits("limitB").abund // 20
                                                        % acc.limits("limitB").id    // 21
                                                        % acc.value("missE")()       // 22
                                                        % acc.value("missI")()       // 23
                                                        % acc.value("missB")()       // 24
                         ).str());
        o.writer->write("\n\n");
    };

    o.writer->open("TransAlign_pooled.stats");
    
    //f(ChrT); TODO
    //f("chr1"); TODO
    

    o.writer->close();
}