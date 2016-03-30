#include "VarQuin/v_coverage.hpp"

using namespace Anaquin;

extern Scripts PlotDensity();

VCoverage::Stats VCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();

    Stats stats;
    
    stats.chrT = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.cID == ChrT ? static_cast<bool>(r.r_var.match(align.l, MatchRule::Contains)) : false;
    });

    /*
     * TODO: Implement statistics for the endogenous. This requires a new standards implementation.
     */

    return stats;
}

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = VCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for the standards
     */

    bo.writer = o.writer;
    bo.file   = "VarCoverage_chrT.bedgraph";

    CoverageTool::bedGraph(stats.chrT, bo, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "VarCoverage_summary.stats";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();

    CoverageTool::summary(stats.chrT, to, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.match(Locus(i, j), MatchRule::Contains);
    });
    
    /*
     * Generating density plot
     */
    
    o.writer->open("VarCoverage_density.R");
    o.writer->write(RWriter::createScript("VarCoverage_chrT.bedgraph", PlotDensity()));
    o.writer->close();
}