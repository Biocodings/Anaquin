#include <catch.hpp>
#include "test.hpp"
#include "RnaQuin/r_express.hpp"

using namespace Anaquin;

TEST_CASE("TExpress_Guided_Invalid")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;
    
    REQUIRE_THROWS(RExpress::analyze("tests/data/guidedInvalid.gtf", o));
    
    Test::clear();
    
    const auto r2 = Test::test("RnaExpression -m data/RnaQuin/MRN027_v001.csv -method gene -ufiles tests/data/guidedInvalid.gtf");
    
    REQUIRE(r2.error == "[ERRO]: Failed to parse \"ABCD\". This is not a number.\n");
    REQUIRE(r2.status == 1);
}

TEST_CASE("TExpress_Replicates")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;
    
    auto r = RExpress::analyze(std::vector<FileName> { "tests/data/A1.gtf", "tests/data/A2.gtf", "tests/data/A3.gtf" }, o);

    REQUIRE(r[0].countSyn   == 170);
    REQUIRE(r[0].countGen   == 0);
    REQUIRE(r[0].dilution() == 1.0);
    REQUIRE(r[1].countSyn   == 3);
    REQUIRE(r[1].countGen   == 0);
    REQUIRE(r[1].dilution() == 1.0);
    REQUIRE(r[2].countSyn   == 170);
    REQUIRE(r[2].countGen   == 0);
    REQUIRE(r[2].dilution() == 1.0);

    REQUIRE(r[0].gData.size() == 0);
    REQUIRE(r[0].genes.size() == 0);
    REQUIRE(r[0].isos.size()  == 136);
    REQUIRE(r[1].gData.size() == 0);
    REQUIRE(r[1].genes.size() == 0);
    REQUIRE(r[1].isos.size()  == 0);
    REQUIRE(r[2].gData.size() == 0);
    REQUIRE(r[2].genes.size() == 0);
    REQUIRE(r[2].isos.size()  == 127);

    REQUIRE(!r[0].isos.count("R2_71_1"));
    REQUIRE(r[0].isos["R2_72_1"].x == Approx(0.0157356262));
    REQUIRE(r[0].isos["R2_72_1"].y == Approx(0.0000099647));
    REQUIRE(r[0].isos["R2_72_2"].x == Approx(0.031471252));
    REQUIRE(r[0].isos["R2_72_2"].y == Approx(0.0000000463));
    REQUIRE(r[0].isos["R2_72_3"].x == Approx(0.062942505));
    REQUIRE(r[0].isos["R2_72_3"].y == Approx(0.0000000047));
    REQUIRE(r[0].isos["R2_72_4"].x == Approx(0.12588501));
    REQUIRE(r[0].isos["R2_72_4"].y == Approx(0.0237112164));

    REQUIRE(!r[1].isos.count("R2_71_1"));
    REQUIRE(!r[1].isos.count("R2_72_1"));
    REQUIRE(!r[1].isos.count("R2_72_2"));
    REQUIRE(!r[1].isos.count("R2_72_3"));
    REQUIRE(!r[1].isos.count("R2_72_4"));
    
    REQUIRE(!r[2].isos.count("R2_71_1"));
    REQUIRE(!r[2].isos.count("R2_72_1"));
    REQUIRE(!r[2].isos.count("R2_72_2"));
    REQUIRE(!r[2].isos.count("R2_72_3"));
    REQUIRE(!r[2].isos.count("R2_72_4"));
    
    REQUIRE(r[0].limit.id == "R2_72_1");
    REQUIRE(r[0].limit.abund == Approx(0.015735626));
    REQUIRE(r[1].limit.id.empty());
    REQUIRE(isnan(r[1].limit.abund));
    REQUIRE(r[2].limit.id == "R1_91_1");
    REQUIRE(r[2].limit.abund == Approx(0.014305115));
}

TEST_CASE("TExpress_Guided_Equal")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;
    
    auto r = RExpress::analyze("tests/data/guidedEqual.gtf", o);
    
    REQUIRE(r.countSyn == 8);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);
    
    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.genes.size() == 0);
    REQUIRE(r.isos.size()  == 8);

    REQUIRE(r.isos["R2_71_1"].x == Approx(0.114440918));
    REQUIRE(r.isos["R2_71_1"].y == Approx(0.05));
    REQUIRE(r.isos["R2_71_2"].x == Approx(3.662109375));
    REQUIRE(r.isos["R2_71_2"].y == Approx(0.05));
    
    REQUIRE(r.isos["R2_72_1"].x == Approx(0.0157356262));
    REQUIRE(r.isos["R2_72_1"].y == Approx(0.05));
    REQUIRE(r.isos["R2_72_2"].x == Approx(0.031471252));
    REQUIRE(r.isos["R2_72_2"].y == Approx(0.05));
    REQUIRE(r.isos["R2_72_3"].x == Approx(0.062942505));
    REQUIRE(r.isos["R2_72_3"].y == Approx(0.05));
    REQUIRE(r.isos["R2_72_4"].x == Approx(0.12588501));
    REQUIRE(r.isos["R2_72_4"].y == Approx(0.05));

    REQUIRE(r.isos["R2_73_1"].x == Approx(0.057220459));
    REQUIRE(r.isos["R2_73_1"].y == Approx(0.05));
    REQUIRE(r.isos["R2_73_2"].x == Approx(1.831054688));
    REQUIRE(r.isos["R2_73_2"].y == Approx(0.05));
    
    REQUIRE(r.limit.id == "R2_72_1");
    REQUIRE(r.limit.abund == Approx(0.015735626));

    const auto ls = r.genes.linear();
    
    REQUIRE(isnan(ls.r));
    REQUIRE(isnan(ls.R2));
    
    Test::clear();
    
    const auto r2 = Test::test("RnaExpression -m data/RnaQuin/MRN027_v001.csv -method gene -ufiles tests/data/guidedEqual.gtf");
    
    REQUIRE(r2.error == "");
    REQUIRE(r2.status == 0);
}

TEST_CASE("TExpress_Guided_Head")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Gene;
    
    auto r = RExpress::analyze("tests/data/guidedHead.gtf", o);
    
    REQUIRE(r.countSyn == 3);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);
    
    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.genes.size() == 3);

    REQUIRE(r.genes["R2_71"].x == Approx(3.776550293));
    REQUIRE(r.genes["R2_71"].y == Approx(0.2522223943));
    REQUIRE(r.genes["R2_72"].x == Approx(0.2360343933));
    REQUIRE(r.genes["R2_72"].y == Approx(0.4431061925));
    REQUIRE(r.genes["R2_73"].x == Approx(1.8882751465));
    REQUIRE(r.genes["R2_73"].y == Approx(0.3761108));
    
    REQUIRE(r.limit.id == "R2_72");
    REQUIRE(r.limit.abund == Approx(0.2360343933));
    
    const auto ls = r.genes.linear();

    REQUIRE(ls.r  == Approx(-0.8688405903));
    REQUIRE(ls.R2 == Approx(0.7548839714));
}

TEST_CASE("TExpress_Denovo_Genes")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Gene;
    
    auto r1 = RExpress::analyze("tests/data/denovo.gtf", o);
    
    REQUIRE(r1.countSyn == 0);
    REQUIRE(r1.countGen == 293);
    REQUIRE(r1.dilution() == 0.0);
    
    REQUIRE(r1.gData.size() == 293);
    REQUIRE(r1.genes.size() == 0);
    REQUIRE(r1.isos.size()  == 0);

    REQUIRE(r1.limit.id.empty());
    REQUIRE(isnan(r1.limit.abund));
    
    Test::clear();
    
    const auto r2 = Test::test("RnaExpression -m data/RnaQuin/MRN027_v001.csv -method gene -ufiles tests/data/denovo.gtf");
    
    REQUIRE(r2.error == "[ERRO]: Failed to find anything on the in-silico chromosome: tests/data/denovo.gtf\n");
    REQUIRE(r2.status == 1);
}

TEST_CASE("TExpress_Denovo_Isoforms")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;

    auto r1 = RExpress::analyze("tests/data/denovo.gtf", o);
    
    REQUIRE(r1.countSyn == 0);
    REQUIRE(r1.countGen == 293);
    REQUIRE(r1.dilution() == 0.0);

    REQUIRE(r1.gData.size() == 293);
    REQUIRE(r1.genes.size() == 0);
    REQUIRE(r1.isos.size()  == 0);
    
    REQUIRE(r1.limit.id.empty());
    REQUIRE(isnan(r1.limit.abund));
    
    Test::clear();
    
    const auto r2 = Test::test("RnaExpression -m data/RnaQuin/MRN027_v001.csv -method isoform -ufiles tests/data/denovo.gtf");
    
    REQUIRE(r2.error == "[ERRO]: Failed to find anything on the in-silico chromosome: tests/data/denovo.gtf\n");
    REQUIRE(r2.status == 1);
}

TEST_CASE("TExpress_Guided_Genes")
{
    Test::transA();
 
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Gene;

    auto r = RExpress::analyze("tests/data/guided.gtf", o);

    REQUIRE(r.countSyn == 74);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);

    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.genes.size() == 74);

    REQUIRE(r.genes["R2_73"].x == Approx(1.8882751465));
    REQUIRE(r.genes["R2_73"].y == Approx(0.3761108));

    const auto ls = r.genes.linear();

    REQUIRE(ls.r  == Approx(0.9520079819));
    REQUIRE(ls.R2 == Approx(0.9063191975));
}

TEST_CASE("TExpress_Guided_Isoforms")
{
    Test::transA();
    
    auto o = RExpress::Options();
    
    o.format = RExpress::Format::GTF;
    o.metrs  = RExpress::Metrics::Isoform;
    
    auto r = RExpress::analyze("tests/data/guided.gtf", o);
    
    REQUIRE(r.countSyn == 214);
    REQUIRE(r.countGen == 0);
    REQUIRE(r.dilution() == 1.0);
    
    REQUIRE(r.gData.size() == 0);
    REQUIRE(r.isos.size() == 149);

    REQUIRE(r.isos["R2_73_1"].x == Approx(0.057220459));
    REQUIRE(r.isos["R2_73_1"].y == Approx(0.0000587338));
    REQUIRE(r.isos["R2_73_2"].x == Approx(1.831054688));
    REQUIRE(r.isos["R2_73_2"].y == Approx(0.3760520666));

    const auto ls = r.isos.linear();
    
    REQUIRE(ls.r  == Approx(0.8052543963));
    REQUIRE(ls.R2 == Approx(0.6484346427));
}