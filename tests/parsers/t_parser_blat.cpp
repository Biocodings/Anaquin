#include <catch.hpp>
#include "parsers/parser_blat.hpp"

using namespace Anaquin;

TEST_CASE("ParserBlat_1")
{
    REQUIRE(ParserBlat::isBlat(Reader("tests/data/Contigs.psl")));
    REQUIRE(!ParserBlat::isBlat(Reader("tests/data/transcripts.gtf")));
}
