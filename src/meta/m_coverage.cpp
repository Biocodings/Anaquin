#include "meta/m_coverage.hpp"

using namespace Anaquin;

MCoverage::Stats MCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    const auto &r = Standard::instance().r_meta;
    
    return CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return r.match(align.id) && r.match(align.id)->l.contains(align.l);
    });
}

void MCoverage::report(const FileName &file, const MCoverage::Options &o)
{
    const auto &r    = Standard::instance().r_meta;
    const auto stats = MCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */
    
    bo.writer = o.writer;
    bo.file   = "MetaCoverage_community.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics for each specie
     */
    
    const auto inters = r.intervals();
    
    for (const auto &i : inters.map())
    {
        CoverageTool::CoverageReportOptions to;
        
        to.writer  = o.writer;
        to.summary = "MetaCoverage_" + i.first + ".stats";
        to.refs    = r.hist().size();
        to.length  = r.size();
        to.id      = i.first;
        
        CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
        {
            // Filter to the regions in the standards
            return r.match(Locus(i, j), MatchRule::Contains);
        });
    }

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "MetaCoverage_summary.stats";
    to.refs     = r.hist().size();
    to.length   = r.size();

/*
 TODO: Need to give a fake name
    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_meta.match(Locus(i, j), MatchRule::Contains);
    });
*/
}