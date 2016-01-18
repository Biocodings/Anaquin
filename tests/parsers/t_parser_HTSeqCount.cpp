#include <catch.hpp>
#include "parsers/parser_HTSeqCount.hpp"

using namespace Anaquin;

TEST_CASE("ParserHTSeqCount_Multiple_Test")
{
    const auto x = "ENSG00000000003.14\t0\n"
                   "ENSG00000000005.5\t3\n"
                   "ENSG00000000419.12\t6";

    const auto y = "ENSG00000000003.14\t1\n"
                   "ENSG00000000005.5\t4\n"
                   "ENSG00000000419.12\t7";

    const auto z = "ENSG00000000003.14\t2\n"
                   "ENSG00000000005.5\t5\n"
                   "ENSG00000000419.12\t8";

    std::vector<ParserHTSeqCount::Samples> samps;
    
    ParserHTSeqCount::parse(std::vector<Reader> { Reader(x, DataMode::String),
                                                  Reader(y, DataMode::String),
                                                  Reader(z, DataMode::String) },
                            [&](const ParserHTSeqCount::Samples &s, const ParserProgress &)
                            {
                                samps.push_back(s);
                            });
    
    REQUIRE(samps.size() == 3);

    REQUIRE(samps[0].id == "ENSG00000000003.14");
    REQUIRE(samps[0].counts[0] == 0);
    REQUIRE(samps[0].counts[1] == 1);
    REQUIRE(samps[0].counts[2] == 2);
    
    REQUIRE(samps[1].id == "ENSG00000000005.5");
    REQUIRE(samps[1].counts[0] == 3);
    REQUIRE(samps[1].counts[1] == 4);
    REQUIRE(samps[1].counts[2] == 5);

    REQUIRE(samps[2].id == "ENSG00000000419.12");
    REQUIRE(samps[2].counts[0] == 6);
    REQUIRE(samps[2].counts[1] == 7);
    REQUIRE(samps[2].counts[2] == 8);
}

TEST_CASE("ParserHTSeqCount_Single_Test")
{
    const auto str = "ENSG00000000003.14\t0\n"
                     "ENSG00000000005.5\t0\n"
                     "ENSG00000000419.12\t71";
    
    std::vector<ParserHTSeqCount::Sample> x;
    
    ParserHTSeqCount::parse(Reader(str, DataMode::String), [&](const ParserHTSeqCount::Sample &r, const ParserProgress &)
    {
        x.push_back(r);
    });
    
    REQUIRE(x.size() == 3);
    
    REQUIRE(x[0].id    == "ENSG00000000003.14");
    REQUIRE(x[0].count == 0);
    REQUIRE(x[1].id    == "ENSG00000000005.5");
    REQUIRE(x[1].count == 0);
    REQUIRE(x[2].id    == "ENSG00000000419.12");
    REQUIRE(x[2].count == 71);
}