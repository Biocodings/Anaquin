#include "meta/m_coverage.hpp"

using namespace Anaquin;

MCoverage::Stats MCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    const auto &r = Standard::instance().r_meta;
    
    return CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return r.match(align.cID) && r.match(align.cID)->l.contains(align.l);
    });
}

void MCoverage::report(const FileName &file, const MCoverage::Options &o)
{
    const auto &r    = Standard::instance().r_meta;
    const auto stats = MCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the community
     */
    
    bo.writer = o.writer;
    bo.file   = "MetaCoverage_community.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating detailed statistics for each MetaQuin
     */

    const std::string format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.writer->open("MetaCoverage_quins.csv");
    o.writer->write((boost::format(format) % "ID"
                                           % "Length"
                                           % "Min"
                                           % "Max"
                                           % "Mean"
                                           % "P25"
                                           % "P50"
                                           % "P75"
                                           % "Covered").str());

    const auto inters = r.intervals();
    
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
    
    /*
     * Generating summary statistics
     */
/*
    const auto summary = "Summary for input: %1%\n\n"
                         "   Experiment: %2%\n"
                         "   Synthetic: %3%\n\n"
                         "   ************ References ************\n\n"
                         "   Sequins: %4%\n"
                         "   Bases:   %5%\n\n"
                         "   ************ Statistics ************\n\n"
                         "   Minimum: %6%\n"
                         "   Maximum: %7%\n"
                         "   Mean:    %8%\n"
                         "   Covered  %9%\n";
*/
    const auto overallStats = inters.stats();

    
    
    
    
    
    
    
    
    
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