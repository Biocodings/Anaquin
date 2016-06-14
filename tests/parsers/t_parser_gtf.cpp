#include <catch.hpp>
#include "parsers/parser_gtf.hpp"

using namespace Anaquin;

TEST_CASE("ParserGTF_FPKM")
{
    std::map<TransID, FPKM> g2f;
    std::map<TransID, FPKM> t2f;
    
    ParserGTF::parse(Reader("tests/data/transcripts.gtf"), [&](
                const ParserGTF::Data &x, const std::string &, const ParserProgress &)
    {
        if (x.type == RNAFeature::Transcript)
        {
            t2f[x.tID] = x.fpkm;
        }
    });
    
    REQUIRE(t2f.size() == 165);
    REQUIRE(t2f["R2_76_1"] == Approx(3.54458));
    REQUIRE(t2f["R2_76_2"] == Approx(1.84317));
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