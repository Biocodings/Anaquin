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
            stats.hist[align.mapped ? align.cID : "NA"]++;
        }

        stats.update(align);

        if (align.mapped && f(align, info.p))
        {
            if (!stats.inters.find(align.cID))
            {
                // Add a new interva for the chromosome
                stats.inters.add(Interval(align.cID, Locus(0, info.length-1)));
            }

            stats.inters.find(align.cID)->add(align.l);
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

    const auto iStats  = inter->stats(f);
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Fraction of reads mapped to the synthetic and experimental chromosomes\n"
                         "   ***\n\n"
                         "   Unmapped:   %2%\n"
                         "   Synthetic:  %3%\n"
                         "   Experiment: %4%\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   Synthetic: %5% sequins\n"
                         "   Synthetic: %6% bases\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Experiment)\n"
                         "   ***\n\n"
                         "   Experiment: %7%\n\n"
                         "   ****************************************************\n"
                         "   ***                                              ***\n"
                         "   ***    Statistics for the coverge (Synthetic)    ***\n"
                         "   ***                                              ***\n"
                         "   ****************************************************\n\n"
                         "   Minimum: %8%\n"
                         "   Maximum: %9%\n"
                         "   Mean:    %10%\n"
                         "   25th:    %11%\n"
                         "   50th:    %12%\n"
                         "   75th:    %13%\n";

    o.writer->open(o.summary);
    o.writer->write((boost::format(summary) % stats.src
                                            % stats.unmapped
                                            % stats.n_chrT
                                            % stats.n_endo
                                            % o.refs
                                            % o.length
                                            % "NA"
                                            % iStats.min
                                            % iStats.max
                                            % iStats.mean
                                            % iStats.p25
                                            % iStats.p50
                                            % iStats.p75).str());
    o.writer->close();
}