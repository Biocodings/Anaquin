#include <catch.hpp>
#include "data/reader.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

TEST_CASE("ParserQuast_Test")
{
    std::vector<ParserQuast::GenomeData> x;
    
    ParserQuast::parseGenomeInfo(Reader("tests/data/genome_info.txt"), [&](const ParserQuast::GenomeData &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x[0].id      == "MG_29");
    REQUIRE(x[0].total   == 2974);
    REQUIRE(x[0].covered == 0);
    REQUIRE(x[1].id      == "M3_G");
    REQUIRE(x[1].total   == 1824);
    REQUIRE(x[1].covered == 1793);
    
    REQUIRE(x.size() == 63);
}

TEST_CASE("ParserQuast_Invalid")
{
    REQUIRE_THROWS(ParserQuast::parseGenomeInfo(Reader("tests/data/DESeq2.csv"),
                    [&](const ParserQuast::GenomeData &d, const ParserProgress &) {}));
}