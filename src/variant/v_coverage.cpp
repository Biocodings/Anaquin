#include "variant/v_coverage.hpp"

using namespace Anaquin;

VCoverage::Stats VCoverage::stats(const FileName &file, const Options &o)
{
    o.analyze(file);
    
    const auto &r = Standard::instance();

    // Reuse a generic tool for coverage
    const auto stats = CoverageTool::stats(file, [&](const Alignment &align, const ParserProgress &)
    {
        if (align.id != r.id)
        {
            return false;
        }

        return static_cast<bool>(r.r_var.findGeno(align.l));
    });

    return stats;
}

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    const auto &r    = Standard::instance();
    const auto stats = VCoverage::stats(file, o);

    CoverageTool::CoverageReportOptions to;

    to.writer   = o.writer;
    to.summary  = "VarCoverage_summary.stats";
    to.bedGraph = "VarCoverage_summary.bedgraph";
    to.refs     = r.r_var.hist().size();
    to.length   = r.r_var.size();

    // Reuse a generic tool for coverage
    CoverageTool::report(stats, to, [&](const ChromoID &id, Base i, Base j, Coverage)
    {
        if (id != r.id)
        {
            return false;
        }

        return static_cast<bool>(r.r_var.findGeno(Locus(i, j)));
    });
}