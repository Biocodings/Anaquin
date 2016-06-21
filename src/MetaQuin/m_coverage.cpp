#include "MetaQuin/m_coverage.hpp"

using namespace Anaquin;

MCoverage::Stats MCoverage::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    const auto &r = Standard::instance().r_meta;

    return CoverageTool::stats_(file, r.hist(), [&](const Alignment &align, const ParserProgress &)
    {
        if (align.cID == ChrT)
        {
            return r.match(align.l, MatchRule::Contains);
        }
        
        return (const SequinData *) nullptr;
        //return r.match(align.cID) && r.match(align.cID)->l.contains(align.l);
    });
}

void MCoverage::report(const FileName &file, const MCoverage::Options &o)
{
    const auto &r    = Standard::instance().r_meta;
    const auto stats = MCoverage::analyze(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * 1. Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "MetaCoverage_summary.stats";
    to.refs     = r.hist().size();
    to.length   = r.size();
    
//    CoverageTool::summary(stats, to, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        // Filter to the regions in the standards
//        return r.match(Locus(i, j), MatchRule::Contains);
//    });
    
    /*
     * 2. Generating bedgraph for the community
     */
    
    bo.writer = o.writer;
    bo.file   = "MetaCoverage_coverage.bedgraph";

//    CoverageTool::bedGraph(stats, bo, [&](const ChrID &id, Base i, Base j, Coverage)
//    {
//        // Filter to the regions in the standards
//        return r.match(Locus(i, j), MatchRule::Contains);
//    });

    /*
     * 3. Generating detailed statistics for each MetaQuin
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->open("MetaCoverage_sequins.csv");
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Min"
                                           % "Max"
                                           % "Mean"
                                           % "P25"
                                           % "P50"
                                           % "P75"
                                           % "Covered").str());

    const auto inters = r.inters();
    
    for (const auto &i : inters.data())
    {
        const auto &stats = i.second.stats();
        
        o.writer->write((boost::format(format) % i.first
                                               % stats.length
                                               % stats.min
                                               % stats.max
                                               % stats.mean
                                               % stats.p25
                                               % stats.p50
                                               % stats.p75
                                               % stats.covered()).str());
    }
    
    o.writer->close();
}