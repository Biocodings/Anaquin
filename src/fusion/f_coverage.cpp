#include "fusion/f_coverage.hpp"

using namespace Anaquin;

FCoverage::Stats FCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();
    Stats stats;
    
    stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.id == Standard::chrT ? static_cast<bool>(r.r_fus.match(align.l, MatchRule::Contains)) : false;
    });

    return stats;
}

void FCoverage::report(const FileName &file, const FCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = FCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */
    
    bo.writer = o.writer;
    bo.file   = "FusionCoverage_chrT.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_fus.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "FusionCoverage_summary.stats";
    to.refs     = r.r_fus.hist().size();
    to.length   = r.r_fus.size();

    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_fus.match(Locus(i, j), MatchRule::Contains);
    });
}