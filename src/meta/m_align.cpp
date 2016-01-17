#include "meta/m_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

MAlign::Stats MAlign::analyze(const FileName &file, const Options &o)
{
    MAlign::Stats stats;

    const auto &r = Standard::instance().r_meta;

    // Intervals for the reference community
    stats.inters = r.intervals();

    o.analyze(file);

    stats.bp.h = stats.sp.h = Standard::instance().r_meta.hist();

    // False positives for each specie
    std::map<GenomeID, Base> fps;

    for (const auto &i : stats.inters.data())
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
            /*
             * Calculating at the sequin level
             */
            
            if (classify(stats.sp.m, align, [&](const Alignment &)
            {
                return match->l().contains(align.l);
            }))
            {
                stats.sp.h.at(match->id())++;
            };
            
            /*
             * Calculating at the base level. Delay the confusion matrix because we'd have to
             * consider overlapping.
             */

            fps.at(align.id) += match->map(align.l);
        }
    });

    o.info("Calculating references");

    // Use the distribution to compute the references
    stats.sp.inferRefFromHist();

    /*
     * Calculating metrics for all references at the base level.
     */
    
    for (const auto &i : stats.inters.data())
    {
        auto &m  = stats.base.at(i.first);
        auto &in = i.second;

        // Update the overall performance
        stats.bp.m.fp() += fps.at(i.first);

        in.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth)
            {
                // Update the sequin performance
                m.tp() += j - i;
                
                // Update the overall performance
                stats.bp.m.tp() += j - i;

                stats.bp.h.at(id)++;
            }
        });

        m.nr() = in.l().length();
        m.nq() = m.tp() + fps.at(i.first);
        
        assert(m.nr() >= m.tp());
        
        stats.bp.m.nr() += in.l().length();
        stats.bp.m.nq()  = stats.bp.m.tp() + stats.bp.m.fp();
    }

    o.info("Calculating detection limit");
    
    stats.bp.limit = r.limit(stats.bp.h);
    stats.sp.limit = r.limit(stats.sp.h);

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
                         "   Unmapped:    %2% alignments\n"
                         "   Community:   %3% alignments\n"
                         "   Synthetic:   %4% alignments\n\n"
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
                                            % stats.n_endo
                                            % stats.n_chrT
                                            % stats.inters.size()                       // 5
                                            % stats.sp.m.sn()
                                            % stats.sp.m.ac()
                                            % stats.sp.limit.id
                                            % stats.sp.limit.abund
                                            % stats.bp.m.sn()    // 10
                                            % stats.bp.m.ac()
                                            % stats.bp.limit.id
                                            % stats.bp.limit.abund
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

    for (const auto &i : stats.inters.data())
    {
        const auto &b = stats.base.at(i.first);
        const auto ss = i.second.stats();

        o.writer->write((boost::format(format) % i.first
                                               % stats.sp.h.at(i.first)
                                               % stats.bp.h.at(i.first)
                                               % ss.covered()
                                               % (isnan(b.ac()) ? "-" : std::to_string(b.ac()))
                        ).str());
    }

    o.writer->close();
}