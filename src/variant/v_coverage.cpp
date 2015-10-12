#include "variant/v_coverage.hpp"

using namespace Anaquin;

void VCoverage::report(const FileName &file, const VCoverage::Options &o)
{
    o.analyze(file);
    
    CoverageAnalyzer::report(file, "VarCoverage_summary.bedgraph", o);
}