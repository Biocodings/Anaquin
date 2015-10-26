#include "ladder/l_coverage.hpp"

using namespace Anaquin;

LCoverage::Stats LCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();

    return CoverageTool::stats(file,  [&](const Alignment &align, const ParserProgress &)
    {
        return align.id == r.id ? static_cast<bool>(r.r_lad.match(align.l, MatchRule::Contains)) : false;
    });
}

void LCoverage::report(const FileName &file, const LCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = LCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */
    
    bo.writer = o.writer;
    bo.file   = "LadderCoverage_chrT.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_lad.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "LadderCoverage_summary.stats";
    to.refs     = r.r_lad.hist().size();
    to.length   = r.r_lad.size();

    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_lad.match(Locus(i, j), MatchRule::Contains);
    });
}