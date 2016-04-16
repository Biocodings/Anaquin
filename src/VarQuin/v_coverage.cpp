#include "VarQuin/v_coverage.hpp"

using namespace Anaquin;

// Defined in resources.cpp
extern Scripts PlotVDensity();

VCoverage::Stats VCoverage::analyze(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();

    return CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.cID == ChrT ? static_cast<bool>(r.r_var.match(align.l, MatchRule::Contains)) : false;
    });
}

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = VCoverage::analyze(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating summary statistics
     */
    
    o.info("Generating VarCoverage_summary.stats");

    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "VarCoverage_summary.stats";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();
    
    CoverageTool::summary(stats, to, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating bedgraph for the standards
     */

    bo.writer = o.writer;
    bo.file   = "VarCoverage_density.bedgraph";

    CoverageTool::bedGraph(stats, bo, [&](const ChrID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.match(Locus(i, j), MatchRule::Contains);
    });

    /*
     * Generating density plot
     */
    
    o.writer->open("VarCoverage_density.R");
    o.writer->write(RWriter::createScript("VarCoverage_chrT.bedgraph", PlotVDensity()));
    o.writer->close();
}