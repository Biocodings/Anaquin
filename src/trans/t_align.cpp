#include <iostream>
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
    stats.pb.h = stats.histE = stats.histI = r.geneHist();

    stats.eInters = r.exonInters();
    stats.iInters = r.intronInters();

    // For each gene...
    for (const auto &i : stats.pb.h)
    {
        stats.sb[i.first];
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
    
    if (inters.contains(align.l, cMatches))
    {
        for (auto &i : cMatches)
        {
            contains.at(i->id())++;
        }
    }
    else
    {
        // Maybe we can find an exon that overlaps the alignment?
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
    
    for (const auto &i : stats.pb.h)
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

    assert(stats.pb.m.tp() == 0 && stats.pb.m.fp() == 0);

    /*
     * 1. Calculating alignment statistics.
     */
    
    o.info("Calculating alignment statistics");

    auto aligns = [](std::map<GeneID, TAlign::MergedConfusion> &gene,
                     TAlign::MergedConfusion &over,
                     Hist &h,
                     Counts unknowns,
                     const TAlign::BinCounts &contains,
                     const TAlign::BinCounts &overlaps,
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
        auto &m  = stats.sb.at(i.second.gID);
        auto &in = i.second;
        
        const auto &gID = i.second.gID;

        // Update the FP at the gene level
        m.fp() = lFPS.at(gID) + rFPS.at(gID);
        
        // Update the FP at the overall level
        stats.pb.m.fp() += m.fp();
        
        Base covered = 0;
        
        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                covered += j - i;
                
                // Update the overall performance
                stats.pb.m.tp() += j - i;
                            
                // Update the distribution
                stats.pb.h.at(gID)++;
            }
        });

        m.tp() += covered;
        m.nr() += in.l().length();
        m.nq()  = m.tp() + m.fp();
        
        assert(m.nr() >= m.tp());
        
        stats.pb.m.nr() += in.l().length();
        stats.pb.m.nq()  = stats.pb.m.tp() + stats.pb.m.fp();
    }

    o.info("Base (TP): " + std::to_string(stats.pb.m.tp()));
    o.info("Base (FP): " + std::to_string(stats.pb.m.fp()));

    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    /*
     * Calculaing the hard detection limit. The limit is estimated with the frequency distribution.
     */
    
    const auto &r = Standard::instance().r_trans;

    stats.limitE   = r.limitGene(stats.histE);
    stats.limitI   = r.limitGene(stats.histI);
    stats.pb.limit = r.limitGene(stats.pb.h);
    
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

void TAlign::report(const FileName &file, const Options &o)
{
    const auto &r = Standard::instance().r_trans;
    const auto stats = TAlign::analyze(file, o);
    
    o.logInfo((boost::format("Base: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pb.m.nr()
                                          % stats.pb.m.nq()
                                          % stats.pb.m.tp()
                                          % stats.pb.m.fp()
                                          % stats.pb.m.fn()
                                          % stats.pb.m.sn()
                                          % stats.pb.m.ac()).str());

    /*
     * Write out summary statistics
     */

    {
        const auto summary = "Summary for dataset: %1%\n\n"
                             "   Unmapped:   %2% (%24%) reads\n"
                             "   Experiment: %3% (%25%) reads\n"
                             "   Synthetic:  %4% reads\n\n"
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
        
        o.writer->open("TransAlign_summary.stats");
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
                                                % stats.sn(Stats::AlignMetrics::Exon)
                                                % stats.precise(Stats::AlignMetrics::Exon)
                                                % stats.limitE.abund
                                                % stats.limitE.id
                                                % stats.sn(Stats::AlignMetrics::Intron)
                                                % stats.precise(Stats::AlignMetrics::Intron)
                                                % stats.limitI.abund
                                                % stats.limitI.id
                                                % stats.sn(Stats::AlignMetrics::Base)
                                                % stats.precise(Stats::AlignMetrics::Base)
                                                % stats.pb.limit.abund
                                                % stats.pb.limit.id
                                                % stats.dilution()
                                                % stats.exp()
                                                % stats.chrT()
                         ).str());
        o.writer->close();
    }

    /*
     * Generating detailed statistics for each sequin
     */
    
    {
        o.writer->open("TransAlign_quins.stats");
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
        
        for (const auto &i : stats.pb.h)
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

            const auto &mb = stats.sb.at(i.first);
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
}