#include <catch.hpp>
#include "tools/samtools.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

TEST_CASE("HT_Test2")
{
    auto f = sam_open("tests/data/test1.bam", "r");
    auto t = bam_init1();
    auto h = sam_hdr_read(f);
    
    int i = 0;
    
    std::vector<std::string> r1;
    
    while (sam_read1(f, h, t) >= 0)
    {
        r1.push_back(bam2cigar(t));

        if (i++ == 24)
        {
            // No need to test all alignments in the file
            break;
        }
    }

    REQUIRE(r1[0]  == "61M");
    REQUIRE(r1[1]  == "61M");
    REQUIRE(r1[2]  == "40M");
    REQUIRE(r1[24] == "44M3D80M");
}
