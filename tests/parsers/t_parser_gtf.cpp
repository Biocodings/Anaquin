#include <catch.hpp>
#include "parsers/parser_gtf.hpp"

using namespace Anaquin;

#ifdef INTERNAL_TESTING

TEST_CASE("ParserGTF_Gencode")
{
    std::vector<Feature> fs;
    
    ParserGTF::parse(Reader("tests/data/GeneCodeV23Annotation.gtf"), [&](const Feature &f, const std::string &, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs.size() == 2228240);
}

#endif