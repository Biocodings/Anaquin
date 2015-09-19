#include "variant/v_align.hpp"
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

VAlign::Stats VAlign::report(const std::string &file, const Options &o)
{
    VAlign::Stats stats;

    const auto &r = Standard::instance().r_var;

    o.info("Parsing alignment file");

    ParserSAM::parse(file, [&](const Alignment &align, const ParserProgress &p)
    {
        if (!align.i && (p.i % 1000000) == 0)
        {
            o.wait(std::to_string(p.i));
        }
        
        if (!align.i)
        {
            if      (!align.mapped)                       { stats.unmapped++; }
            else if (align.id != Standard::instance().id) { stats.n_hg38++;   }
            else                                          { stats.n_chrT++;   }
        }
        
        if (!align.mapped || align.id != Standard::instance().id)
        {
            return;
        }

        /*
         * Collect statistics for the alignment
         */
        
        const SequinData * match;

        if (classify(stats.p.m, align, [&](const Alignment &)
        {
            return (match = r.match(align.l, Contains)) ? Positive : Negative;
        }))
        {
            stats.h.at(match->id)++;
        }
    });

    sums(stats.h, stats.p.m.nr);

    o.info("Calculating limit of sensitivity");

    // Calculate for the sensitivity
    stats.p.s = r.limit(stats.h);

    o.info("Generating statistics");

    /*
     * Write out summary statistics
     */
    
    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Unmapped:    %2% reads\n"
                         "   Genome:      %3% reads\n"
                         "   Synthetic:   %4% reads\n\n"
                         "   Reference:   %5% genes\n"
                         "   Sensitivity: %6%\n"
                         "   Specificity: %7%\n\n"
                         "   Dilution:    %8%\n"
    ;

    o.writer->open("VarAlign_summary.stats");
    o.writer->write((boost::format(summary) % file
                                            % stats.unmapped
                                            % stats.n_hg38
                                            % stats.n_chrT
                                            % (r.countRefGenes() + r.countVarGens())
                                            % stats.p.m.sn()
                                            % stats.p.m.sp()
                                            % stats.dilution()).str());
    o.writer->close();
    
    /*
     * Write out sequin statistics
     */
    
    o.writer->open("VarAlign_quins.stats");
    o.writer->write((boost::format("Summary for dataset: %1%\n") % file).str());
    
    const auto format = "%1%\t%2%";
    o.writer->write((boost::format(format) % "id" % "counts (reads)").str());
    
    for (const auto &i : stats.h)
    {
        o.writer->write((boost::format(format) % i.first % stats.h.at(i.first)).str());
    }
    
    o.writer->close();

	return stats;
}
