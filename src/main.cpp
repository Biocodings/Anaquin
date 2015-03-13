#include "align.hpp"

int main(int argc, char * const argv[])
{
    
    AlignAnalyzer::analyze("accepted_hits.sam");
   
    
    struct AlignAnalyzer
    {
        /*
         * Analyze sequence alignments for sequins for a given SAM file.
         */
        
        static void analyze(const std::string &file);

    
    
    return 0;
}