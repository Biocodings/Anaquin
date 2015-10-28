#include "meta/m_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

template <typename T, typename Classifer> bool classifyBase(Confusion &m, const T &t, Classifer c)
{
    const auto data = c(t);

    if (data)
    {
        const auto overs = data->l().overlap(t.l);

        /*
         * Anything overlaps is a true-positive, while everything else is a false negative.
         */

        m.nq   += t.l.length();
        m.tp() += overs;
        m.fp() += (m.nq - m.tp());
    }

    return data;
}

MAlign::Stats MAlign::analyze(const FileName &file, const Options &o)
{
    MAlign::Stats stats;

    const auto &r = Standard::instance().r_meta;

    stats.inters = r.intervals();
    o.analyze(file);

    stats.p[PerfLevel::BasePerf].h   = Standard::instance().r_meta.hist();
    stats.p[PerfLevel::SequinPerf].h = Standard::instance().r_meta.hist();

    // False positives for each specie
    std::map<GenomeID, Base> fps;

    for (const auto &i : stats.inters.map())
    {
        fps[i.first] = 0;
        stats.base[i.first];
    }

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        REPORT_STATUS();

        // Does this alignment belongs to one of the synthetic species?
        Interval * const match = stats.inters.find(align.id);

        stats.update(align, [&](const Alignment &)
        {
            return match;
        });

        if (!align.mapped)
        {
            return;
        }

        if (match)
        {
            auto &sp = stats.p.at(PerfLevel::SequinPerf);

            /*
             * Calculating at the sequin level
             */
            
            if (classify(sp.m, align, [&](const Alignment &)
            {
                return match->l().contains(align.l);
            }))
            {
                sp.h.at(match->id())++;
            };
            
            /*
             * Calculating at the base level. Delay the confusion matrix because we'd have to
             * consider overlapping.
             */

            fps.at(align.id) += match->add_(align.l);
        }
    });

    o.info("Calculating references");
    
    auto sp = &(stats.p.at(PerfLevel::SequinPerf));
    auto bp = &(stats.p.at(PerfLevel::BasePerf));

    sums(sp->h, sp->m.nr);
    
    /*
     * Metrics at the base level for all reference genomes
     */
    
    bp->m.tp() = bp->m.fp() = bp->m.nq = 0;
    
    for (const auto &i : stats.inters.map())
    {
        auto &m  = stats.base.at(i.first);
        auto &in = i.second;

        // Update the overall performance
        bp->m.fp() += fps.at(i.first);
        
        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                m.tp() += j - i;
                
                // Update the overall performance
                bp->m.tp() += j - i;

                bp->h.at(id)++;
            }
        });

        m.nr = in.l().length();
        m.nq = m.tp() + fps.at(i.first);
        
        assert(m.nr >= m.tp());
        
        bp->m.nr += in.l().length();
        bp->m.nq  = bp->m.tp() + bp->m.fp();
    }

    o.info("Calculating detection limit");
    
    bp->s = r.limit(bp->h);
    sp->s = r.limit(sp->h);

    return stats;
}

void MAlign::report(const FileName &file, const Options &o)
{
    const auto stats = MAlign::analyze(file, o);
    
    /*
     * Generating summary statistics
     */

    o.info("Generating summary statistics");

    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:    %2% reads\n"
                         "   Community:   %3% reads\n"
                         "   Synthetic:   %4% reads\n\n"
                         "   Reference:   %5% sequins\n\n"
                         "   ************ Sequin Level ************\n\n"
                         "   Sensitivity:     %6%\n"
                         "   Accuracy:        %7%\n"
                         "   Detection Limit: %8% (%9%)\n\n"
                         "   ************ Base Level ************\n\n"
                         "   Sensitivity:     %10%\n"
                         "   Accuracy:        %11%\n\n"
                         "   Detection limit: %12% (%13%)\n\n"
                         "   Dilution: %14%\n";

    o.writer->open("MetaAlign_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.unmapped
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % stats.inters.size()                       // 5
                                            % stats.p.at(PerfLevel::SequinPerf).m.sn()
                                            % stats.p.at(PerfLevel::SequinPerf).m.sp()
                                            % stats.p.at(PerfLevel::SequinPerf).s.id
                                            % stats.p.at(PerfLevel::SequinPerf).s.abund
                                            % stats.p.at(PerfLevel::BasePerf).m.sn()    // 10
                                            % stats.p.at(PerfLevel::BasePerf).m.sp()
                                            % stats.p.at(PerfLevel::BasePerf).s.id
                                            % stats.p.at(PerfLevel::BasePerf).s.abund
                                            % stats.dilution()).str());
    o.writer->close();

    /*
     * Generating detailed statistics for each sequin
     */
    
    o.writer->open("MetaAlign_quins.stats");
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());

    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    o.writer->write((boost::format(format) % "ID"
                                           % "Counts (Sequin)"
                                           % "Counts (Base)"
                                           % "Sensitivity"
                                           % "Accuracy"
                     ).str());

    for (const auto &i : stats.inters.map())
    {
        const auto &b  = stats.base.at(i.first);
        const auto &bh = stats.p.at(PerfLevel::BasePerf).h;
        const auto &sh = stats.p.at(PerfLevel::SequinPerf).h;
        
        const auto a = i.second;
        const auto ss  = i.second.stats();

        o.writer->write((boost::format(format) % i.first
                                               % sh.at(i.first)
                                               % bh.at(i.first)
                                               % ss.covered()
                                               % (isnan(b.sp()) ? "-" : std::to_string(b.sp()))
                        ).str());
    }

    o.writer->close();
}
