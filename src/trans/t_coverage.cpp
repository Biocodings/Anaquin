#include "trans/t_coverage.hpp"

using namespace Anaquin;

TCoverage::Stats TCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();
    Stats stats;
    
    stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.id == Standard::chrT ? static_cast<bool>(r.r_trans.match(align.l, MatchRule::Contains)) : false;
    });

    return stats;
}

void TCoverage::report(const FileName &file, const TCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = TCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */
    
    bo.writer = o.writer;
    bo.file   = "TransCoverage_chrT.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_trans.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "TransCoverage_summary.stats";
    to.refs     = r.r_trans.hist().size();
    to.length   = r.r_trans.size();

    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_trans.match(Locus(i, j), MatchRule::Contains);
    });
}