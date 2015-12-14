#include <ss/stats.hpp>
#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

typedef std::map<GeneID, Base> FPStats;

struct ParseImpl
{
    TAlign::Stats *stats;
    
    FPStats  *lFPS;
    FPStats  *rFPS;
    
    Interval *base;
};

// Internal implementation
typedef std::function<void (const ParseImpl &)> Functor;

static TAlign::Stats init()
{
    const auto &r = Standard::instance().r_trans;
    
    TAlign::Stats stats;

    // Initalize the distributions
    stats.overB.h = stats.histE = stats.histI = r.geneHist();

    stats.eInters = r.exonInters();
    stats.iInters = r.intronInters();

    // For each gene...
    for (const auto &i : stats.overB.h)
    {
        stats.geneB[i.first];
        stats.geneE[i.first];
        stats.geneI[i.first];
    }
    
    // For each exon bin...
    for (const auto &i : stats.eInters.data())
    {
        stats.eContains[i.first];
        stats.eOverlaps[i.first];
        stats.exonToGene[i.second.id()] = i.second.gID;
    }
    
    // For each intron bin...
    for (const auto &i : stats.iInters.data())
    {
        stats.iContains[i.first];
        stats.iOverlaps[i.first];
        stats.intronToGene[i.second.id()] = i.second.gID;
    }
    
    assert(!stats.eContains.empty() && !stats.eOverlaps.empty());
    assert(!stats.iContains.empty() && !stats.iOverlaps.empty());
    
    return stats;
}

template <typename T> const T * matchT(const Alignment &align,
                                       Intervals<T> &inters,
                                       std::map<std::string, Counts> &contains,
                                       std::map<std::string, Counts> &overlaps,
                                       FPStats *lFPS = nullptr,
                                       FPStats *rFPS = nullptr)
{
    std::vector<T *> oMatches, cMatches;

    /*
     * It's quite likely there'll be more than a match. Note that it's not possible to distinguish the individuals due to
     * alternative splicing. Thus, we simply increment for all the possible matches. Consequently, it's not possible to detect
     * anything at the isoform level.
     */

    if (inters.contains(align.l, cMatches))
    {
        for (auto &i : cMatches)
        {
            contains.at(i->id())++;
        }
    }
    else
    {
        if (inters.contains(align.l, oMatches))
        {
            for (auto &i : cMatches)
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

TAlign::Stats calculate(const TAlign::Options &o, Functor f)
{
    TAlign::Stats stats = init();
    
    FPStats lFPS, rFPS;
    
    for (const auto &i : stats.overB.h)
    {
        lFPS[i.first];
        rFPS[i.first];
    }
    
    // This is needed to track the FP at the base level
    Interval base("Base", Locus(0, 44566700));
    
    ParseImpl impl;
    
    impl.lFPS  = &lFPS;
    impl.rFPS  = &rFPS;
    impl.base  = &base;
    impl.stats = &stats;

    // It's the caller job to handle the parsing
    f(impl);

    assert(stats.overB.m.tp() == 0 && stats.overB.m.fp() == 0);

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

    aligns(stats.geneE,
           stats.overE,
           stats.histE,
           stats.unknowns.size(),
           stats.eContains,
           stats.eOverlaps,
           stats.exonToGene);

    aligns(stats.geneI,
           stats.overI,
           stats.histI,
           0,
           stats.iContains,
           stats.iOverlaps,
           stats.intronToGene);

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
    
    genes(stats.geneE, stats.overE, stats.eContains, stats.eOverlaps, stats.exonToGene);
    genes(stats.geneI, stats.overI, stats.iContains, stats.iOverlaps, stats.intronToGene);
    
    /*
     * 3. Calculating metrics at the base level.
     */

    o.info("Calculating base statistics");

    for (const auto &i : stats.eInters.data())
    {
        auto &m  = stats.geneB.at(i.second.gID);
        auto &in = i.second;
        
        const auto &gID = i.second.gID;

        // Update the FP at the gene level
        m.fp() = lFPS.at(gID) + rFPS.at(gID);
        
        // Update the FP at the overall level
        stats.overB.m.fp() += m.fp();
        
        Base covered = 0;
        
        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                covered += j - i;
                
                // Update the overall performance
                stats.overB.m.tp() += j - i;
                            
                // Update the distribution
                stats.overB.h.at(gID)++;
            }
        });

        m.tp() += covered;
        m.nr() += in.l().length();
        m.nq()  = m.tp() + m.fp();
        
        assert(m.nr() >= m.tp());
        
        stats.overB.m.nr() += in.l().length();
        stats.overB.m.nq()  = stats.overB.m.tp() + stats.overB.m.fp();
    }

    o.info("Base (TP): " + std::to_string(stats.overB.m.tp()));
    o.info("Base (FP): " + std::to_string(stats.overB.m.fp()));

    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    const auto &r = Standard::instance().r_trans;

    stats.limitE = r.limitGene(stats.histE);
    stats.limitI = r.limitGene(stats.histI);
    stats.overB.limit = r.limitGene(stats.overB.h);
    
    /*
     * Calculating for missing statistics
     */
    
    {
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
        missing(stats.missE, stats.eContains);

        // An intron is missing if no alignment aligns to it
        missing(stats.missI, stats.iContains);
        
        /*
         * A gene is considered missing if not all exons have alignment aligned to it
         */
        
        for (const auto &gene : stats.histE)
        {
            bool missing = false;
            
            for (const auto &bin : stats.eContains)
            {
                if (gene.first == stats.exonToGene.at(bin.first))
                {
                    if ((missing = (bin.second == 0)))
                    {
                        break;
                    }
                }
            }
            
            if (missing)
            {
                stats.missG.insert(Missing(gene.first));
            }
        }
    }
    
    return stats;
}

static void update(const ParseImpl &impl, const Alignment &align, const ParserSAM::AlignmentInfo &info, const TAlign::Options &o)
{
    REPORT_STATUS();
    
    impl.stats->update(align);
    
    if (!align.mapped || align.id != Standard::chrT)
    {
        return;
    }
    
    const Interval *match = nullptr;
    
    if (!align.spliced)
    {
        if ((match = matchT(align,
                            impl.stats->eInters,
                            impl.stats->eContains,
                            impl.stats->eOverlaps,
                            impl.lFPS,
                            impl.rFPS)))
        {
            o.logInfo((boost::format("Exon (match): %1% %2%") % align.id % match->id()).str());
        }
    }
    else
    {
        if ((match = matchT(align,
                            impl.stats->iInters,
                            impl.stats->iContains,
                            impl.stats->iOverlaps)))
        {
            o.logInfo((boost::format("Intron (match): %1% %2%") % align.id % match->id()).str());
        }
    }

    if (!match)
    {
        impl.stats->unknowns.push_back(UnknownAlignment(align.qName, align.l));
        
        // We can't simply add it to statistics because we'll need to account for overlapping
        impl.base->map(align.l);
    }
}

TAlign::Stats TAlign::analyze(const std::vector<Alignment> &aligns, const Options &o)
{
    return calculate(o, [&](const ParseImpl &impl)
    {
        ParserSAM::AlignmentInfo info;
        
        for (const auto &align : aligns)
        {
            update(impl, align, info, o);
        }
    });
}

TAlign::Stats TAlign::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    return calculate(o, [&](const ParseImpl &impl)
    {
        ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
        {
            update(impl, align, info, o);
        });
    });
}

std::vector<TAlign::Stats> TAlign::analyze(const std::vector<FileName> &files, const Options &o)
{
    std::vector<TAlign::Stats> stats;
    
    for (const auto &file : files)
    {
        stats.push_back(analyze(file, o));
    }
    
    return stats;
}

std::vector<TAlign::Stats> TAlign::analyze(const std::vector<std::vector<Alignment>> &aligns, const Options &o)
{
    std::vector<TAlign::Stats> stats;
    
    for (const auto &align : aligns)
    {
        stats.push_back(analyze(align, o));
    }
    
    return stats;
}

static std::string summary()
{
    return "Summary for dataset: %1%\n\n"
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
           "   Detection:   %24% (%25%)\n";
}

// Write summary statistics for a single replicate
static void writeSummary(const TAlign::Stats &stats, const FileName &file, const TAlign::Options &o)
{
    const auto &r = Standard::instance().r_trans;

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:   %2% reads\n"
                         "   Experiment: %3% (%24%%%) reads\n"
                         "   Synthetic:  %4% (%25%%%) reads\n\n"
                         "   Reference:  %5% exons\n"
                         "   Reference:  %6% introns\n"
                         "   Reference:  %7% bases\n\n"
                         "   Query:      %8% exons\n"
                         "   Query:      %9% introns\n"
                         "   Query:      %10% bases\n\n"
                         "   Dilution:   %23%\n\n"
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
                         "   Sensitivity: %11%\n"
                         "   Specificity: %12%\n"
                         "   Detection:   %13% (%14%)\n\n"
                         "   -------------------- Intron level --------------------\n\n"
                         "   Sensitivity: %15%\n"
                         "   Specificity: %16%\n"
                         "   Detection:   %17% (%18%)\n\n"
                         "   -------------------- Base level --------------------\n\n"
                         "   Sensitivity: %19%\n"
                         "   Specificity: %20%\n"
                         "   Detection:   %21% (%22%)\n";
    
    o.writer->open(file);
    o.writer->write((boost::format(summary) % file
                                            % stats.unmapped
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % r.countSortedExons()
                                            % r.countSortedIntrons()
                                            % r.exonBase()
                                            % stats.qExons()
                                            % stats.qIntrons()
                                            % stats.qBases()
                                            % stats.sn(TAlign::Stats::AlignMetrics::AlignExon)
                                            % stats.pc(TAlign::Stats::AlignMetrics::AlignExon)
                                            % stats.limitE.abund
                                            % stats.limitE.id
                                            % stats.sn(TAlign::Stats::AlignMetrics::AlignIntron)
                                            % stats.pc(TAlign::Stats::AlignMetrics::AlignIntron)
                                            % stats.limitI.abund
                                            % stats.limitI.id
                                            % stats.sn(TAlign::Stats::AlignMetrics::AlignBase)
                                            % stats.pc(TAlign::Stats::AlignMetrics::AlignBase)
                                            % stats.overB.limit.abund
                                            % stats.overB.limit.id
                                            % stats.dilution()
                                            % (100.0 * stats.exp())
                                            % (100.0 * stats.chrT())
                     ).str());
    o.writer->close();
}

// Write sequin statistics for a single replicate
static void writeSequins(const TAlign::Stats &stats, const FileName &file, const TAlign::Options &o)
{
    o.writer->open(file);
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    
    o.writer->write((boost::format(format) % "ID"
                                           % "Covered"
                                           % "Sensitivity (Exon)"
                                           % "Specificity (Exon)"
                                           % "Sensitivity (Intron)"
                                           % "Specificity (Intron)"
                                           % "Sensitivity (Base)"
                                           % "Specificity (Base)").str());
    
    for (const auto &i : stats.overB.h)
    {
        Base length   = 0;
        Base nonZeros = 0;
        
        for (const auto &j : stats.eInters.data())
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
        
        const auto &mb = stats.geneB.at(i.first);
        const auto &me = stats.geneE.at(i.first);
        const auto &mi = stats.geneI.at(i.first);
        
        // Not all sequins have an intron...
        if (mi.lNR)
        {
            o.writer->write((boost::format(format) % i.first
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
            o.writer->write((boost::format(format) % i.first
                                                   % me.sn()
                                                   % me.precise()
                                                   % "--"
                                                   % "--"
                                                   % mb.sn()
                                                   % mb.ac()
                                                   % covered).str());
        }
    }
    
    o.writer->close();
}

template <typename T> class Accumulator
{
    public:
    
        typedef std::string Key;
    
        struct Deviation
        {
            // First moment
            T mean;

            // Standard deviation
            T sd;

            inline std::string operator()() const
            {
                return (boost::format("%1% \u00B1 %2%") % mean % sd).str();
            }
        };
    
        void add(const Key &key, double value)
        {
            _data[key].push_back(value);
        }

        void add(const Key &key, const Limit &l)
        {
            _limits[key].push_back(l);
        }
    
        Deviation value(const std::string &key) const
        {
            Deviation d;
            
            d.sd   = SS::sd(_data.at(key));
            d.mean = SS::mean(_data.at(key));
            
            return d;
        }
    
        const Limit & limits(const Key &key) const
        {
            const Limit *min = nullptr;
        
            for (const auto &limit : _limits.at(key))
            {
                if (!min || limit.abund < min->abund)
                {
                    min = &limit;
                }
            }

            return *min;
        }

    private:
    
        // Used for comparing limit of detection
        std::map<Key, std::vector<Limit>> _limits;
    
        // Used for regular mappings
        std::map<Key, std::vector<double>> _data;
};

void TAlign::report(const std::vector<FileName> &files, const Options &o)
{
    const auto &r = Standard::instance().r_trans;

    const auto stats = TAlign::analyze(files, o);
    
    for (auto i = 0; i < files.size(); i++)
    {
        writeSequins(stats[i], (boost::format("TransAlign_%1%_quins.stats")   % files[i]).str(), o);
        writeSummary(stats[i], (boost::format("TransAlign_%1%_summary.stats") % files[i]).str(), o);
    }

    /*
     * Write pooled summary statistics
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
    
    Accumulator<double> acc;
    
    for (const auto &stat : stats)
    {
        acc.add("Unmapped",    stat.unmapped);
        acc.add("Experiment",  stat.n_expT);
        acc.add("Synthetic",   stat.n_chrT);
        acc.add("QExon",       stat.qExons());
        acc.add("QIntron",     stat.qIntrons());
        acc.add("QBase",       stat.qBases());
        acc.add("Dilution",    stat.n_chrT);
        acc.add("ExonSN",      stat.sn(TAlign::Stats::AlignMetrics::AlignExon));
        acc.add("ExonPC",      stat.pc(TAlign::Stats::AlignMetrics::AlignExon));
        acc.add("IntronSN",    stat.sn(TAlign::Stats::AlignMetrics::AlignIntron));
        acc.add("IntronPC",    stat.pc(TAlign::Stats::AlignMetrics::AlignIntron));
        acc.add("BaseSN",      stat.sn(TAlign::Stats::AlignMetrics::AlignBase));
        acc.add("BasePC",      stat.pc(TAlign::Stats::AlignMetrics::AlignBase));
        acc.add("ExpPercent",  100.0 * stat.exp());
        acc.add("ChrTPercent", 100.0 * stat.chrT());
        acc.add("LimitE",      stat.limitE);
        acc.add("LimitI",      stat.limitI);
        acc.add("LimitB",      stat.overB.limit);
    }
    
    o.writer->open("TransAlign_summary.stats");
    o.writer->write((boost::format(summary()) % concated
                                              % acc.value("Unmapped")()
                                              % acc.value("Experiment")()
                                              % acc.value("ExpPercent")()  // 4
                                              % acc.value("Synthetic")()
                                              % acc.value("ChrTPercent")() // 6
                                              % r.countSortedExons()
                                              % r.countSortedIntrons()
                                              % r.exonBase()
                                              % acc.value("QExon")()
                                              % acc.value("QIntron")()
                                              % acc.value("QBase")()
                                              % acc.value("Dilution")() // 13
                                              % acc.value("ExonSN")()   // 14
                                              % acc.value("ExonPC")()   // 15
                                              % acc.limits("LimitE").abund
                                              % acc.limits("LimitE").id
                                              % acc.value("IntronSN")() // 18
                                              % acc.value("IntronPC")() // 19
                                              % acc.limits("LimitI").abund
                                              % acc.limits("LimitI").id
                                              % acc.value("BaseSN")()   // 22
                                              % acc.value("BasePC")()   // 23
                                              % acc.limits("LimitB").abund
                                              % acc.limits("LimitB").id
                     ).str());
    o.writer->close();
}

void TAlign::report(const FileName &file, const Options &o)
{
    const auto stats = TAlign::analyze(file, o);
    
    writeSummary(stats, "TransAlign_summary.stats", o);
    writeSequins(stats, "TransAlign_quins.stats", o);
}