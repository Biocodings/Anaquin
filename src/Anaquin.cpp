#include "gtest/gtest.h"
#include "AlignerAnalyst.hpp"
#include "AssemblyAnalyst.hpp"
#include <boost/program_options/cmdline.hpp>

int main(int argc, char ** argv)
{
    
    
    
    
//    AlignerAnalyst::analyze("/Users/tedwong/Sources/ABCD/aligned_output/accepted_hits.sam");
//    AlignerAnalyst::analyze("/Users/tedwong/Sources/ABCD/aligned_output/accepted_hits.sam", 1000);
    

    
//    AssemblyAnalyst::analyze("/Users/tedwong/Sources/ABCD/transcripts/transcripts.gtf");
    
  
   // return 0;
    
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}