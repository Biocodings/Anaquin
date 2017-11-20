#include <catch.hpp>
#include "parsers/parser_bam2.hpp"

using namespace Anaquin;

TEST_CASE("ParserBAM2_1")
{
    std::vector<ParserBAM2::Data> r;
    
    ParserBAM2::parse("tests/data/insert.sam", [&](ParserBAM2::Data &x, const ParserBAM2::Info &)
    {
        r.push_back(x);
    });
    
    REQUIRE(r.size() == 1);
    REQUIRE(r[0].cID == "chrT");
    REQUIRE(r[0].l.start == 8288747);
    REQUIRE(r[0].l.end   == 8288769);
}

TEST_CASE("ParserBAM2_2")
{
    std::vector<ParserBAM2::Data> r;

    ParserBAM2::parse("tests/data/deletion.sam", [&](ParserBAM2::Data &x, const ParserBAM2::Info &i)
    {
        r.push_back(x);
    });
    
    REQUIRE(r.size() == 2);

    REQUIRE(r[0].cID == "chrT");
    REQUIRE(r[0].l.start == 7058914);
    REQUIRE(r[0].l.end   == 7059037);
    REQUIRE(r[1].cID == "chrT");
    REQUIRE(r[1].l.start == 7058781);
    REQUIRE(r[1].l.end   == 7058838);
}
