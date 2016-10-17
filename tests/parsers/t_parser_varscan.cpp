#include <catch.hpp>
#include "data/reader.hpp"
#include "parsers/parser_varscan.hpp"

using namespace Anaquin;

TEST_CASE("ParserVarScan_Pile")
{
    std::vector<ParserVarScan::Data> x;
    
    ParserVarScan::parse(Reader("tests/data/varscan.tab"), [&](const ParserVarScan::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x.size() == 1443);
}

TEST_CASE("ParserVarScan_SomaticSNP")
{
    std::vector<ParserVarScan::Data> x;
    
    ParserVarScan::parse(Reader("tests/data/VarQuin.snp"), [&](const ParserVarScan::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x.size() == 139);
}

TEST_CASE("ParserVarScan_SomaticIndel")
{
    std::vector<ParserVarScan::Data> x;
    
    ParserVarScan::parse(Reader("tests/data/VarQuin.indel"), [&](const ParserVarScan::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x.size() == 62);
}
