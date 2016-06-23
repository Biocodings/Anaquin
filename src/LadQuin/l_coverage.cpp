#include "LadQuin/l_coverage.hpp"

using namespace Anaquin;

LCoverage::Stats LCoverage::stats(const FileName &file, const Options &o)
{
    //    const auto &r = Standard::instance().r_lad;

    o.analyze(file);

    throw "Not Implemented";
    
//    return CoverageTool::stats_(file, r.hist(), [&](const Alignment &align, const ParserProgress &)
//    {
//        if (align.cID == ChrT)
//        {
//            return r.match(align.l, MatchRule::Contains);
//        }
//        
//        return (const SequinData *) nullptr;
//    });
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

//    CoverageTool::bedGraph(stats, bo, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        // Filter to the regions in the standards
//        return r.r_lad.match(Locus(i, j), MatchRule::Contains);
//    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "LadderCoverage_summary.stats";
    to.refs     = r.r_lad.hist().size();
    to.length   = r.r_lad.size();

//    CoverageTool::summary(stats, to, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        // Filter to the regions in the standards
//        return r.r_lad.match(Locus(i, j), MatchRule::Contains);
//    });
}