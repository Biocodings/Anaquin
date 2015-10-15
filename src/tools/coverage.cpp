#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, AlignFunctor f)
{
    CoverageTool::Stats stats;

    stats.src = file;
    
    /*
     * Reference: https://github.com/arq5x/bedtools2/blob/master/src/genomeCoverageBed/genomeCoverageBed.cpp
     */

    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i)
        {
            if      (!align.mapped)                       { stats.unmapped++; }
            else if (align.id != Standard::instance().id) { stats.n_expT++;   }
            else                                          { stats.n_chrT++;   }
        }

        // Proceed with the alignment?
        if (f(align, info.p))
        {
            if (!stats.inters.find(align.id))
            {
                stats.inters.add(Interval(align.id, Locus(0, info.size-1)));
            }

            stats.inters.find(align.id)->add(align.l);
        }
    });

    return stats;
}

void CoverageTool::report(const CoverageTool::Stats &stats, const CoverageReportOptions &o, CoverageFunctor f)
{
    o.writer->open(o.bedGraph);

    // Coverage for the synthetic chromosome
    const auto synth = stats.inters.find(Standard::instance().id);

    if (!synth)
    {
        throw std::runtime_error("Failed to find coverage for the synthetic chromosome");
    }

    /*
     * Generating a bedgraph of coverage
     */
    
    synth->bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
    {
        if (depth)
        {
            o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id
                                                                 % i
                                                                 % j
                                                                 % depth).str());
        }
    });

    o.writer->close();

    /*
     * Generating summary statistics
     */

    const auto sstats  = synth->stats(f);
    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Experiment: %2%\n"
                         "   Synthetic: %3%\n\n"
                         "   Reference: %4%\n"
                         "   Reference Bases: %5%\n\n"
                         "   Minimum: %6%\n"
                         "   Maximum: %7%\n\n"
                         "   Minimum: %6%\n"
                         "   Maximum: %7%\n"
                         "   Mean:    %8%\n"
                         "   25th: %9%\n"
                         "   50th: %10%\n"
                         "   75th: %11%\n"
    ;

    o.writer->open(o.summary);
    o.writer->write((boost::format(summary) % stats.src
                                            % stats.n_expT
                                            % stats.n_chrT
                                            % o.refs
                                            % o.length
                                            % sstats.min
                                            % sstats.max
                                            % sstats.mean
                                            % sstats.p25
                                            % sstats.p50
                                            % sstats.p75
                     ).str());
    o.writer->close();
}