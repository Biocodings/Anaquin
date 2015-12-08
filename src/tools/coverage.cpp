#include "tools/coverage.hpp"
#include "parsers/parser_sam.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, AlignFunctor f)
{
    CoverageTool::Stats stats;

    stats.src = file;
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        if (!align.i)
        {
            stats.hist[align.mapped ? align.id : "NA"]++;
        }
        
        stats.update(align);

        if (align.mapped && f(align, info.p))
        {
            if (!stats.inters.find(align.id))
            {
                stats.inters.add(Interval(align.id, Locus(0, info.length-1)));
            }

            stats.inters.find(align.id)->add(align.l);
        }
    });

    return stats;
}

CoverageTool::Stats stats(const std::vector<CoverageTool::Mapping> &maps, const std::map<GenericID, Base> &ref)
{
    CoverageTool::Stats stats;
    
    for (const auto &m : maps)
    {
        if (!stats.inters.find(m.id))
        {
            stats.inters.add(Interval(m.id, Locus(0, ref.at(m.id))));
        }
        
        stats.inters.find(m.id)->add(m.l);
    }

    return stats;
}

void CoverageTool::bedGraph(const Stats &stats, const CoverageBedGraphOptions &o, CoverageFunctor f)
{
    o.writer->open(o.file);

    for (const auto &i : stats.inters.data())
    {
        const auto chr = i.second;

        chr.bedGraph([&](const ChromoID &id, Base i, Base j, Base depth)
        {
            if (depth && f(id, i, j, depth))
            {
                o.writer->write((boost::format("%1%\t%2%\t%3%\t%4%") % id
                                                                     % i
                                                                     % j
                                                                     % depth).str());
            }
        });
    }

    o.writer->close();
}

void CoverageTool::summary(const CoverageTool::Stats &stats, const CoverageReportOptions &o, CoverageFunctor f)
{
    const auto inter = stats.inters.find(o.id);

    if (!inter)
    {
        return;
    }

    /*
     * Generating summary statistics
     */

    const auto sstats  = inter->stats(f);
    const auto summary = "Summary for dataset: %1%\n\n"
                         "   Experiment: %2%\n"
                         "   Synthetic: %3%\n\n"
                         "   Reference: %4%\n"
                         "   Reference Bases: %5%\n\n"
                         "   Minimum: %6%\n"
                         "   Maximum: %7%\n"
                         "   Mean:    %8%\n"
                         "   25th: %9%\n"
                         "   50th: %10%\n"
                         "   75th: %11%\n";

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
                                            % sstats.p75).str());
    o.writer->close();
}