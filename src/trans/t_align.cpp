#include "trans/t_align.hpp"
#include "parsers/parser_sam.hpp"
#include <iostream>
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
    stats.pb.h = stats.pe.h = stats.pi.h = stats.alignExon.h = stats.alignIntron.h = r.geneHist();

    stats.eInters = r.exonInters();
    stats.iInters = r.intronInters();

    // For each gene...
    for (const auto &i : stats.pe.h)
    {
        stats.sb[i.first];
        stats.si[i.first];
        stats.se[i.first];
        stats.detectExons[i.first];
        stats.undetectExons[i.first];
        stats.detectIntrons[i.first];
        stats.undetectIntrons[i.first];
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

    /*
     * Initalize the references
     */
    
    for (auto &i : stats.se)
    {
        for (const auto &j : stats.eInters.data())
        {
            if (i.first == j.second.gID)
            {
                i.second.nr()++;
            }
        }
    }
    
    for (auto &i : stats.si)
    {
        for (const auto &j : stats.iInters.data())
        {
            if (i.first == j.second.gID)
            {
                i.second.nr()++;
            }
        }
    }
    
    return stats;
}

template <typename T> const T * matchT(const Alignment &align,
                                       Intervals<T> &inters,
                                       std::map<std::string, Counts> &contains,
                                       std::map<std::string, Counts> &overlaps,
                                       FPStats *lFPS = nullptr,
                                       FPStats *rFPS = nullptr)
{
    T * oMatch = nullptr;

    // Can we find an exon that contains the alignment?
    T * cMatch = inters.contains(align.l);

    if (cMatch)
    {
        contains.at(cMatch->id())++;
    }
    else
    {
        // Maybe we can find an exon that overlaps the alignment?
        oMatch = inters.overlap(align.l);

        if (oMatch)
        {
            overlaps.at(oMatch->id())++;
        }
    }
    
    T * match = cMatch ? cMatch : oMatch;
    
    if (match)
    {
        Base lp, rp;
        
        // Anything that fails to being mapped is counted as FP
        match->map(align.l, &lp, &rp);

        if (lFPS && rFPS)
        {
            lFPS->at(match->gID) = std::max(lFPS->at(match->gID), lp);
            rFPS->at(match->gID) = std::max(rFPS->at(match->gID), rp);
        }
    }
    
    return cMatch;
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
    
    o.info("Calculating for non-overlapping at the base level");
    
    base.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
    {
        if (depth)
        {
            // Everything here has no mapping and therefore FP
            stats.pb.m.fp() += j - i;
        }
    });
    
    /*
     * 1. Calculating alignment statistics. Those can be used for accuracy at the exon and intron level.
     */
    
    o.info("Calculating alignment statistics");

    auto aligns = [](Counts mapped, Counts unknown, const std::map<ExonID, Counts> &overlaps, Performance &p)
    {
        p.m.tp() = mapped;
        p.m.fp() = unknown;

        for (const auto &i : overlaps)
        {
            p.m.fp() += i.second;
        }
    };

    aligns(stats.eMapped, stats.eUnknown, stats.eOverlaps, stats.alignExon);
    aligns(stats.iMapped, stats.iUnknown, stats.iOverlaps, stats.alignIntron);
    
    /*
     * 2. Calculating metrics at the base level.
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
    
    /*
     * 3. Calculating statistics for sequins (at the gene level due to alternative splicing)
     *
     * Exon:
     *
     *    TP -> for all alignments contained
     *    FP -> for all alignments overlapped
     *    FN -> if TP is zero (the gene is undetected)
     *
     * Intron is similar.
     */
    
    o.info("Calculating statistics for sequins");

    auto sequins = [](std::map<GeneID, Confusion> &m,
                      std::map<std::string, std::string> &mapper,
                      std::map<GeneID, Counts> &detects,
                      std::map<GeneID, Counts> &undetects,
                      const std::map<std::string, Counts> &contains,
                      const std::map<std::string, Counts> &overlaps)
    {
        for (auto &i : m)
        {
            Counts detect = 0;
            
            for (auto &j : contains)
            {
                if (i.first == mapper.at(j.first))
                {
                    /*
                     * Eg: If there're 100 alignments contained, they're all TP.
                     *
                     *    j.second == 100
                     */
                    
                    i.second.tp() += j.second;
                    
                    /*
                     * How to define detection? It's detected if there is at least a single contained
                     * alignment. This can be potentially enhanced to a custom rule.
                     */
                    
                    if (j.second)
                    {
                        detect++;
                        detects.at(i.first)++;
                    }
                }
            }
            
            for (auto &j : overlaps)
            {
                if (i.first == mapper.at(j.first))
                {
                    i.second.fp() += j.second;
                }
            }

            /*
             * How many exons/introns mapped to the gene?
             */
            
            Counts n = 0;
            
            for (const auto &j : mapper)
            {
                if (i.first == j.second)
                {
                    n++;
                }
            }
            
            // Obviously we can't detect more than what we have
            assert(detect <= n);
            
            i.second.fn() = n - detect;
            i.second.nr() = i.second.tp() + i.second.fn();

            undetects.at(i.first) = n - detect;
        }
    };
    
    sequins(stats.se, stats.exonToGene, stats.detectExons, stats.undetectExons, stats.eContains, stats.eOverlaps);
    sequins(stats.si, stats.intronToGene, stats.detectIntrons, stats.undetectIntrons, stats.iContains, stats.iOverlaps);

    /*
     * 4. Calculating overall metrics for exons and introns.
     *
     * Exon:
     *
     *    TP -> for each detected exon
     *    FN -> for each undetected exon
     *
     * This is done for measuring the sensitivity at exon level. For example, if an experiment
     * is only able to detect half of the exons, the sensitivity would be 50%.
     *
     * Unfortunately, there is no concept of FP here. Therefore, FP is undefined. Intron is similar.
     */
    
    o.info("Calculating overall statistics");

    assert(!stats.pe.m.tp() && !stats.pe.m.fp() && !stats.pe.m.fn());
    assert(!stats.pi.m.tp() && !stats.pi.m.fp() && !stats.pi.m.fn());

    // Combinating statistics for exons
    for (const auto &i : stats.se)
    {
        stats.pe.m.tp() += i.second.tp();
        stats.pe.m.fp() += i.second.fp();
        stats.pe.m.fn() += i.second.fn();
        stats.pe.m.nq() += i.second.nq();
        stats.pe.m.nr() += i.second.nr();
    }

    // Combinating statistics for introns
    for (const auto &i : stats.si)
    {
        stats.pi.m.tp() += i.second.tp();
        stats.pi.m.fp() += i.second.fp();
        stats.pi.m.fn() += i.second.fn();
        stats.pi.m.nq() += i.second.nq();
        stats.pi.m.nr() += i.second.nr();
    }

    /*
     * Calculating detection limit
     */
    
    o.info("Calculating detection limit");
    
    /*
     * Calculaing the hard detection limit. The limit is estimated with the frequency distribution.
     */
    
    const auto &r = Standard::instance().r_trans;

    stats.pb.limit          = r.limitGene(stats.pb.h);
    stats.alignExon.limit   = r.limitGene(stats.alignExon.h);
    stats.alignIntron.limit = r.limitGene(stats.alignIntron.h);
    
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
            impl.stats->exonToGene.at(match->id());
            impl.stats->alignExon.h.at(impl.stats->exonToGene.at(match->id()))++;
            impl.stats->eMapped++;
        }
        else
        {
            impl.stats->eUnknown++;
        }
    }
    else
    {
        if ((match = matchT(align,
                            impl.stats->iInters,
                            impl.stats->iContains,
                            impl.stats->iOverlaps)))
        {
            impl.stats->alignIntron.h.at(impl.stats->intronToGene.at(match->id()))++;
            impl.stats->iMapped++;
        }
        else
        {
            impl.stats->iUnknown++;
        }
    }

    if (!match)
    {
        impl.stats->unknowns.push_back(UnknownAlignment(align.qName, align.l));
        
        // We can't simply add it to statistics because we'll need to account for overlapping
        impl.base->map(align.l);
    }
}

TAlign::Stats TAlign::stats(const std::vector<Alignment> &aligns, const Options &o)
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

TAlign::Stats TAlign::stats(const FileName &file, const Options &o)
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
    const auto stats = TAlign::stats(file, o);
    
    o.logInfo((boost::format("Exon: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pe.m.nr()
                                          % stats.pe.m.nq()
                                          % stats.pe.m.tp()
                                          % stats.pe.m.fp()
                                          % stats.pe.m.fn()
                                          % stats.pe.m.sn()
                                          % stats.pe.m.ac()).str());

    o.logInfo((boost::format("Intron: %1% %2% %3% %4% %5% %6% %7%")
                                          % stats.pi.m.nr()
                                          % stats.pi.m.nq()
                                          % stats.pi.m.tp()
                                          % stats.pi.m.fp()
                                          % stats.pi.m.fn()
                                          % stats.pi.m.sn()
                                          % stats.pi.m.ac()).str());

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
                             "   Unmapped:   %2% alignments\n"
                             "   Experiment: %3% alignments\n"
                             "   Synthetic:  %4% alignments\n\n"
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
                             "   Accuarcy:    %12%\n"
                             "   Detection:   %13% (%14%)\n\n"
                             "   -------------------- Intron level --------------------\n\n"
                             "   Sensitivity: %15%\n"
                             "   Accuarcy:    %16%\n"
                             "   Detection:   %17% (%18%)\n\n"
                             "   -------------------- Base level --------------------\n\n"
                             "   Sensitivity: %19%\n"
                             "   Accuarcy:    %20%\n"
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
                                                % stats.ac(Stats::AlignMetrics::Exon)
                                                % stats.alignExon.limit.abund
                                                % stats.alignExon.limit.id
                                                % stats.sn(Stats::AlignMetrics::Intron)
                                                % stats.ac(Stats::AlignMetrics::Intron)
                                                % stats.alignIntron.limit.abund
                                                % stats.alignIntron.limit.id
                                                % stats.sn(Stats::AlignMetrics::Base)
                                                % stats.ac(Stats::AlignMetrics::Base)
                                                % stats.pb.limit.abund
                                                % stats.pb.limit.id
                                                % stats.dilution()).str());
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
                                               % "Accuracy (Exon)"
                                               % "Sensitivity (Intron)"
                                               % "Accuracy (Intron)"
                                               % "Sensitivity (Base)"
                                               % "Accuracy (Base)").str());
        
        for (const auto &i : stats.pe.h)
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

            const auto &me = stats.se.at(i.first);
            const auto &mi = stats.si.at(i.first);
            const auto &mb = stats.sb.at(i.first);

            // Not all sequins have an intron...
            if (mi.nr())
            {
                o.writer->write((boost::format(format) % i.first
                                                       % covered
                                                       % me.sn()
                                                       % me.ac()
                                                       % mi.sn()
                                                       % mi.ac()
                                                       % mb.sn()
                                                       % mb.ac()).str());
            }
            else
            {
                o.writer->write((boost::format(format) % i.first
                                                       % me.sn()
                                                       % me.ac()
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