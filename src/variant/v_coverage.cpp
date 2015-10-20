#include "variant/v_coverage.hpp"

using namespace Anaquin;

VCoverage::Stats VCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();
    Stats stats;
    
    // Calculate coverage for the synthetic chromosome
    stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        return align.id == r.id ? static_cast<bool>(r.r_var.findGeno(align.l)) : false;
    });

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
    bo.chr    = "chrT";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.findGeno(Locus(i, j));
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "VarCoverage_summary.stats";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();

    CoverageTool::summary(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        // Filter to the regions in the standards
        return r.r_var.findGeno(Locus(i, j));
    });
}