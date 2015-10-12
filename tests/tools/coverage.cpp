#include <catch.hpp>
#include "analyzers/coverage.hpp"

using namespace Anaquin;

TEST_CASE("CoverageAnalyzer_Test")
{
    /*
     * bedtools genomecov -ibam intersected.bam -bg
     */
    
    CoverageAnalyzer::Options o;
    
    
    
    CoverageAnalyzer::report("intersected.bam", "abcd.bedgraph", o);
    
    
    
}