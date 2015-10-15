#include "variant/v_coverage.hpp"

using namespace Anaquin;

VCoverage::Stats VCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();

    Stats stats;
    
    // Calculate coverage for chromosome chrT
    stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        if (align.id == "chr21")
        {
            return static_cast<bool>(r.r_var.findInterval("chr21", align.l));
        }
        else if (align.id == r.id)
        {
            return static_cast<bool>(r.r_var.findGeno(align.l));
        }
        
        return false;
    });

    return stats;
}

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = VCoverage::stats(file, o);

    CoverageTool::CoverageBedGraphOptions bo;

    /*
     * Generating bedgraph for synthetic chromosome
     */
    
    bo.writer = o.writer;
    bo.file   = "VarCoverage_chrT.bedgraph";
    bo.chr    = "chrT";

    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        return r.r_var.findGeno(Locus(i, j));
    });

    /*
     * Generating bedgraph for chromosome 21
     */
    
    bo.file = "VarCoverage_chr21.bedgraph";
    bo.chr  = "chr21";
    
    CoverageTool::bedGraph(stats, bo, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        return r.r_var.findInterval("chr21", Locus(i, j));
    });

    /*
     * Generating summary statistics
     */
    
    CoverageTool::CoverageReportOptions to;
    
    to.writer   = o.writer;
    to.summary  = "VarCoverage_summary.stats";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();

    CoverageTool::report(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        return static_cast<bool>(r.r_var.findGeno(Locus(i, j)));
    });
}