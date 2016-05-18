#include "tools/coverage.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

CoverageTool::Stats CoverageTool::stats(const FileName &file, AlignFunctor f)
{
    CoverageTool::Stats stats;

    stats.src  = file;
    
    ParserSAM::parse(file, [&](const Alignment &align, const ParserSAM::AlignmentInfo &info)
    {
        stats.update(align);

        if (align.mapped && f(align, info.p))
        {
            if (!stats.inters.find(align.cID))
            {
                // Add a new interval for the chromosome
                stats.inters.add(Interval(align.cID, Locus(0, info.length-1)));
            }
            
            stats.hist[align.cID]++;
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

        chr.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
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

Scripts CoverageTool::writeCSV(const CoverageTool::Stats &stats, const CoverageReportOptions &o)
{
    const auto format = "%1%\t%2%\n";
    
    std::stringstream ss;
    ss << ((boost::format(format) % "seq" % "count").str());

    for (const auto &seq : stats.hist)
    {
        ss << ((boost::format(format) % seq.first % seq.second).str());
    }

    return ss.str();
}

void CoverageTool::summary(const CoverageTool::Stats &stats,
                           const CoverageReportOptions &o,
                           CoverageFunctor f)
{
    const auto inter = stats.inters.find(o.id);

    if (!inter)
    {
        return;
    }

    const auto iStats  = inter->stats(f);
    const auto summary = "Summary for input: %1%\n\n"
                         "   ***\n"
                         "   *** Proportion of alignments mapped to the synthetic and genome\n"
                         "   ***\n\n"
                         "   Unmapped:  %2%\n"
                         "   Synthetic: %3%\n"
                         "   Genome:    %4%\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Synthetic)\n"
                         "   ***\n\n"
                         "   File: %5%\n\n"
                         "   Synthetic: %6% sequins\n"
                         "   Synthetic: %7% bases\n\n"
                         "   ***\n"
                         "   *** Reference annotation (Genome)\n"
                         "   ***\n\n"
                         "   File: %8%\n\n"
                         "   Genome: %9% intervals\n"
                         "   Genome: %10% bases\n\n"
                         "   ****************************************************\n"
                         "   ***                                              ***\n"
                         "   ***    Statistics for the synthetic chromosome   ***\n"
                         "   ***                                              ***\n"
                         "   ****************************************************\n\n"
                         "   Minimum: %11%\n"
                         "   Maximum: %12%\n"
                         "   Mean:    %13%\n"
                         "   25th:    %14%\n"
                         "   50th:    %15%\n"
                         "   75th:    %16%\n\n"
                         "   ****************************************************\n"
                         "   ***                                              ***\n"
                         "   ***         Statistics for the genome            ***\n"
                         "   ***                                              ***\n"
                         "   ****************************************************\n\n"
                         "   Minimum: %17%\n"
                         "   Maximum: %18%\n"
                         "   Mean:    %19%\n"
                         "   25th:    %20%\n"
                         "   50th:    %21%\n"
                         "   75th:    %22%\n\n";

    o.writer->open(o.summary);
    o.writer->write((boost::format(summary) % stats.src
                                            % stats.unmapped
                                            % stats.n_chrT
                                            % stats.n_geno
                                            % o.rChrT                           // 5
                                            % o.refs                            // 6
                                            % o.length                          // 7
                                            % (o.rGeno.empty() ? "-" : o.rGeno) // 8
                                            % "-"                               // 9
                                            % "-"                               // 10
                                            % iStats.min                        // 11
                                            % iStats.max                        // 12
                                            % iStats.mean                       // 13
                                            % iStats.p25                        // 14
                                            % iStats.p50                        // 15
                                            % iStats.p75                        // 16
                                            % "-"                               // 17
                                            % "-"                               // 18
                                            % "-"                               // 19
                                            % "-"                               // 20
                                            % "-"                               // 21
                                            % "-"                               // 22
                     ).str());
    o.writer->close();
}