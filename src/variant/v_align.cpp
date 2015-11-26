#include "variant/v_align.hpp"
#include "variant/v_sample.hpp"
#include "parsers/parser_sam.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

VAlign::Stats VAlign::report(const FileName &file, const Options &o)
{
    VAlign::Stats stats;

    const auto &r = Standard::instance().r_var;

    o.info("Parsing alignment file");

    Intervals<> ii;
    
    for (auto &i : r.data())
    {
        if (boost::algorithm::ends_with(i.first, "_V"))
        {
            continue;
        }
        
        const auto id = i.first.substr(0, i.first.size() - 2);
        
        Interval in(id, i.second.l);
        std::cout << id << " " << i.second.l.start << " " << i.second.l.end << std::endl;
        ii.add(in);
    }

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i && !(info.p.i % 1000000))
        {
            o.wait(std::to_string(info.p.i));
        }
        
        stats.update(align);
        
        if (!align.mapped || align.id != Standard::chrT)
        {
            return;
        }

        /*
         * Collect statistics for the alignment
         */

        const VarRef::GenotypeData * match;

        if (classify(stats.p.m, align, [&](const Alignment &)
        {
            return (match = r.findGeno(align.l, Contains)) ? Positive : Negative;
        }))
        {
            const auto l = ii.find(match->id);
            Locus t = Locus(align.l.start - l->l().start, align.l.end - l->l().end);
            ii.find(match->id)->add(t);
 
            stats.h.at(match->id)++;
        }
    });

    sums(stats.h, stats.p.m.nr);

    o.info("Calculating limit of sensitivity");

    // Calculate for the sensitivity
    stats.p.s = r.limitGeno(stats.h);

    o.logInfo((boost::format("Performance: %1% %2% %3% %4% %5% %6% %7%")
                                    % stats.p.m.nr
                                    % stats.p.m.nq
                                    % stats.p.m.tp()
                                    % stats.p.m.fp()
                                    % stats.p.m.fn()
                                    % stats.p.m.sn()
                                    % stats.p.m.sp()).str());
    o.info("Generating summary statistics");

    double covered = 0;
    
    for (auto i : ii.map())
    {
        covered += i.second.stats().covered();
        std::cout << i.second.stats().covered() << std::endl;
    }
    
    covered = covered / ii.map().size();
    
    /*
     * Write out summary statistics
     */
    
    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:    %2% reads\n"
                         "   Experiment:  %3% reads\n"
                         "   Synthetic:   %4% reads\n\n"
                         "   Reference:   %5% genes\n"
                         "   Sensitivity: %6%\n"
                         "   Accuracy:    %7%\n\n"
                         "   Base Covered:    %8%\n\n"
                         "   Detection limit: %9% (%10%)\n\n"
                         "   Dilution:    %11%\n";

    o.writer->open("VarAlign_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.unmapped
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % (r.countRefGenes() + r.countVarGens())
                                            % stats.p.m.sn()
                                            % stats.p.m.sp()
                                            % covered
                                            % stats.p.s.abund
                                            % stats.p.s.id
                                            % stats.dilution()).str());
    o.writer->close();
    
    /*
     * Generating detailed statistics for each sequin
     */
    
    o.writer->open("VarAlign_quins.stats");
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "ID" % "Counts (reads)").str());
    
    for (const auto &i : stats.h)
    {
        o.writer->write((boost::format(format) % i.first % stats.h.at(i.first)).str());
    }
    
    o.writer->close();

	return stats;
}
