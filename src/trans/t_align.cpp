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
    stats.pb.h = stats.pe.h = stats.pi.h = r.geneHist();

    for (const auto &i : stats.pe.h)
    {
        const auto gID = i.first;
        
        stats.sb[gID];
        stats.si[gID];
        stats.se[gID];
        stats.detectExons[gID];
        stats.undetectExons[gID];
        stats.detectIntrons[gID];
        stats.undetectIntrons[gID];
    }
    
    stats.eInters = r.exonInters();
    stats.iInters = r.intronInters();

    // For each exon bin...
    for (const auto &i : stats.eInters.data())
    {
        stats.eContains[i.first];
        stats.eOverlaps[i.first];
        stats.exonToGene[i.second.id()] = i.second.gID;
    }
    
    // For each intron bin..
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
    T * match = nullptr;

    // Can we find an exon that contains the alignment?
    match = inters.contains(align.l);

    if (match)
    {
        contains.at(match->id())++;
    }
    else
    {
        // Maybe we can find an exon that overlaps the alignment?
        match = inters.overlap(align.l);

        if (match)
        {
            overlaps.at(match->id())++;
        }
    }
    
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
    
    return match;
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
    
    o.info("Counting references");
    
    stats.pe.inferRefFromHist();
    stats.pi.inferRefFromHist();
    
    /*
     * Calculating statistics for all sequins
     */

    o.info("Calculating statistics for all sequins");

    assert(!stats.pe.m.tp());
    
    /*
     * Calculating overall metrics for the exon and intron level.
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
    
    auto overall = [&](const std::map<ExonID, Counts> &contains, Confusion &m)
    {
        for (const auto &i : contains)
        {
            if (i.second)
            {
                // It's TP because it's contained
                m.tp()++;
            }
            else
            {
                // It's FP because it is overlapped or outside
                m.fn()++;
            }
        }
    };

    overall(stats.eContains, stats.pe.m);
    overall(stats.iContains, stats.pi.m);
    
    stats.pe.m.fp() = stats.eUnknown; // TODO: add overlaps
    stats.pi.m.fp() = stats.iUnknown;

    /*
     * Calculating alignment statistics. They can be used for accuracy at the exon and intron level.
     */
    
    auto aligns = [](Counts mapped, Counts unknown, const std::map<ExonID, Counts> &overlaps, Confusion &m)
    {
        m.tp() = mapped;
        m.fp() = unknown;
        
        for (const auto &i : overlaps)
        {
            m.fp() += i.second;
        }
    };

    aligns(stats.eMapped, stats.eUnknown, stats.eOverlaps, stats.ae);
    aligns(stats.iMapped, stats.iUnknown, stats.iOverlaps, stats.ai);
    
    /*
     * Calculating metrics at the base level.
     */

    o.info("Calculating statistics for each sequin");

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
     * Exon:
     *
     *    TP -> for all alignments contained
     *    FP -> for all alignments overlapped
     *    FN -> if TP is zero (the gene is undetected)
     *
     * Intron is similar.
     */
    
    auto indiv = [](std::map<GeneID, Confusion> &m,
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
                    i.second.tp() += j.second;
                    
                    if (j.second)
                    {
                        detect++;
                        
                        //
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
    
    indiv(stats.se, stats.exonToGene, stats.detectExons, stats.undetectExons, stats.eContains, stats.eOverlaps);
    indiv(stats.si, stats.intronToGene, stats.detectIntrons, stats.undetectIntrons, stats.iContains, stats.iOverlaps);

    o.info("Calculating detection limit");
    
    const auto &r = Standard::instance().r_trans;

    stats.pe.hl = r.limitGene(stats.pe.h);
    stats.pi.hl = r.limitGene(stats.pi.h);
    stats.pb.hl = r.limitGene(stats.pb.h);
    
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
                             "   Unmapped:   %2% reads\n"
                             "   Experiment: %3% reads\n"
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
                             "   *** Exon level is defined by the performance per exon. An alignment that\n"
                             "   *** is not mapped entirely within an exon is considered as a FP. The\n"
                             "   *** intron level is similar.\n"
                             "   ***\n"
                             "   *** Base level is defined by the performance per nucleotide. A partial\n"
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
                                                % stats.pe.m.nq()
                                                % stats.pi.m.nq()
                                                % stats.pb.m.nq()
                                                % stats.pe.m.sn()
                                                % stats.pe.m.ac()
                                                % stats.pe.hl.abund
                                                % stats.pe.hl.id
                                                % stats.pi.m.sn()
                                                % stats.pi.m.ac()
                                                % stats.pi.hl.abund
                                                % stats.pi.hl.id
                                                % stats.pb.m.sn()
                                                % stats.pb.m.ac()
                                                % stats.pb.hl.abund
                                                % stats.pb.hl.id
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
                                               % "Sensitivity (Exon)"
                                               % "Accuracy (Exon)"
                                               % "Sensitivity (Intron)"
                                               % "Accuracy (Intron)"
                                               % "Sensitivity (Base)"
                                               % "Accuracy (Base)"
                                               % "Covered").str());
        
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

            o.writer->write((boost::format(format) % i.first
                                                   % stats.se.at(i.first).sn()
                                                   % stats.se.at(i.first).ac()
                                                   % stats.si.at(i.first).sn()
                                                   % stats.si.at(i.first).ac()
                                                   % stats.sb.at(i.first).sn()
                                                   % stats.sb.at(i.first).ac()
                                                   % covered).str());
        }
        
        o.writer->close();
    }
}