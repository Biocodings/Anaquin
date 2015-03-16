#include "gtest/gtest.h"
#include "AlignerAnalyst.hpp"
#include "AssemblyAnalyst.hpp"

int main(int argc, char ** argv)
{
    
//    AlignerAnalyst::analyze("/Users/user1/Sources/ABCD/aligned_output/accepted_hits.sam");
//    AlignerAnalyst::analyze("/Users/user1/Sources/ABCD/aligned_output/accepted_hits.sam", 1000);
    

    
//    AssemblyAnalyst::analyze("/Users/user1/Sources/ABCD/transcripts/transcripts.gtf");
    
  
   // return 0;
    
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}