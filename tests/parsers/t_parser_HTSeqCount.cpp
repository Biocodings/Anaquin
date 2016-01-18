#include <catch.hpp>
#include "parsers/parser_HTSeqCount.hpp"

using namespace Anaquin;

TEST_CASE("ParserHTSeqCount_Test")
{
    const auto str = "ENSG00000000003.14\t0\n"
                     "ENSG00000000005.5\t0\n"
                     "ENSG00000000419.12\t71";
    
    std::vector<ParserHTSeqCount::CountRow> x;
    
    ParserHTSeqCount::parse(Reader(str, DataMode::String), [&](const ParserHTSeqCount::CountRow &r, const ParserProgress &)
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