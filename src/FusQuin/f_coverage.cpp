#include "FusQuin/f_coverage.hpp"

using namespace Anaquin;

FCoverage::Stats FCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance().r_fus;
    
    return CoverageTool::stats_(file, r.hist(), [&](const Alignment &align, const ParserProgress &)
    {
        if (align.cID == ChrT)
        {
            return r.match(align.l, MatchRule::Contains);
        }
        
        return (const SequinData *) nullptr;
    });
}

void FCoverage::report(const FileName &file, const FCoverage::Options &o)
{
    const auto &r    = Standard::instance().r_fus;
    const auto stats = FCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    o.info("Generating statistics");
    
    /*
     * 1. Generating summary statistics
     */
    
    o.info("Generating FusionCoverage_summary.stats");
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "FusionCoverage_summary.stats";
    to.refs     = r.hist().size();
    to.length   = r.size();
    
    CoverageTool::summary(stats, to, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.match(Locus(i, j), MatchRule::Contains);
    });
    
    /*
     * 2. Generating bedgraph for the standards
     */
    
    o.info("Generating FusionCoverage_coverage.bedgraph");
    
    bo.writer = o.writer;
    bo.file   = "FusionCoverage_coverage.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.match(Locus(i, j), MatchRule::Contains);
    });
}