#include <catch.hpp>
#include "tools/coverage.hpp"

using namespace Anaquin;

TEST_CASE("CoverageTool_Test")
{
    /*
     * bedtools genomecov -ibam intersected.bam -bg
     */
    
    
    
    CoverageTool::report("intersected.bam", "abcd.bedgraph", [&](const Alignment &align, const ParserProgress &p)
    {
        return true;
    });
    
    
    
    
  //  const auto stats = CoverageTool::analyze("intersected.bam", [&](const Alignment &align, const ParserProgress /&p)
    //{
      //  return true;
    //});

    
}