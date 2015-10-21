#include "meta/m_coverage.hpp"

using namespace Anaquin;

MCoverage::Stats MCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();
    Stats stats;
    
    stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.id == r.id ? static_cast<bool>(r.r_meta.match(align.l, MatchRule::Contains)) : false;
    });

    return stats;
}

void MCoverage::report(const FileName &file, const MCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = MCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */
    
    bo.writer = o.writer;
    bo.file   = "MetaCoverage_chrT.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_meta.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "MetaCoverage_summary.stats";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();

    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_meta.match(Locus(i, j), MatchRule::Contains);
    });
}