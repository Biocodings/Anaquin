#include <catch.hpp>
#include "tools/samtools.hpp"
#include "parsers/parser_sam.hpp"

#include <iostream>
using namespace Anaquin;

TEST_CASE("HT_Reverse")
{
    ParserSAM::parse("data/test/VarQuin/test1.bam", [&](ParserSAM::Data &x, const ParserSAM::Info &i)
    {
        if (i.p.i == 0)
        {
            reverse(x, i);

            const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";
            const auto str = (boost::format(format) % x.name
                                                    % x.flag
                                                    % x.cID
                                                    % x.l.start
                                                    % x.mapq
                                                    % x.cigar
                                                    % x.rnext
                                                    % x.pnext
                                                    % x.tlen
                                                    % x.seq
                                                    % x.qual);
            std::cout << str << std::endl;
        }
    });
}

TEST_CASE("HT_Test2")
{
    auto f = sam_open("data/test/VarQuin/test1.bam", "r");
    auto t = bam_init1();
    auto h = sam_hdr_read(f);
    
    int i = 0;
    
    std::vector<std::string> r1;
    std::vector<std::string> r2;
    
    while (sam_read1(f, h, t) >= 0)
    {
        r1.push_back(bam2cigar(t));
        r2.push_back(bam2rcigar(t));
        
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

    REQUIRE(r2[0]  == "61M");
    REQUIRE(r2[1]  == "61M");
    REQUIRE(r2[2]  == "40M");
    REQUIRE(r2[24] == "80M3D44M");
}