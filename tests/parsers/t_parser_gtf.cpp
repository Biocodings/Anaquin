#include <catch.hpp>
#include "parsers/parser_gtf2.hpp"

using namespace Anaquin;

TEST_CASE("ParserGTF_1")
{
    std::map<TransID, ParserGTF2::Data> t2d;
    
    ParserGTF2::parse(Reader("tests/data/transcripts.gtf"), [&](const ParserGTF2::Data &x, const ParserProgress &)
    {
        if (x.type == RNAFeature::Transcript)
        {
            REQUIRE(!x.cID.empty());
            REQUIRE(!isnan(x.fpkm));
            
            t2d[x.tID] = x;
            
            REQUIRE(!t2d[x.tID].cID.empty());
            REQUIRE(!isnan(t2d[x.tID].fpkm));
        }
    });
    
    REQUIRE(t2d.size() == 165);
    REQUIRE(t2d.count("R1_12_1"));
    REQUIRE(t2d.count("R1_63_1"));

    REQUIRE(t2d["R1_12_1"].cID == "chrT");
    REQUIRE(t2d["R1_12_1"].gID == "R1_12");
    REQUIRE(t2d["R1_12_1"].tID == "R1_12_1");
    REQUIRE(t2d["R1_12_1"].l.end == 10233186);
    REQUIRE(t2d["R1_12_1"].l.start == 10218062);
    REQUIRE(t2d["R1_12_1"].fpkm == 2.6592866492);
    
    REQUIRE(t2d["R1_63_1"].cID == "chrT");
    REQUIRE(t2d["R1_63_1"].gID == "R1_63");
    REQUIRE(t2d["R1_63_1"].tID == "R1_63_1");
    REQUIRE(t2d["R1_63_1"].l.end == 9015138);
    REQUIRE(t2d["R1_63_1"].l.start == 8945052);
    REQUIRE(t2d["R1_63_1"].fpkm == 1783.9793649586);
}

#ifdef INTERNAL_TESTING

TEST_CASE("ParserGTF_Gencode")
{
    std::vector<Feature> fs;
    
    ParserGTF::parse(Reader("tests/data/GeneCodeV23Annotation.gtf"), [&](const ParserGTF::Data &f, const std::string &, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs.size() == 2228240);
}

#endif