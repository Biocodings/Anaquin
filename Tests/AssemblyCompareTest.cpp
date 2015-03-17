#include "gtest/gtest.h"
#include "AssemblyAnalyst.hpp"

using namespace std;

TEST(Generated, AssemblyCompareTest)
{
    const auto r = AssemblyAnalyst::analyze("/Users/user1/Sources/ABCD/transcripts/transcripts.gtf");
    
    std::cout << r.exon.tp << std::endl;
    std::cout << r.exon.fp << std::endl;

}