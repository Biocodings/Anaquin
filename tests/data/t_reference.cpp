#include <catch.hpp>
#include "unit/test.hpp"
#include "data/standard.hpp"

using namespace Anaquin;

TEST_CASE("VarQuin_Fold")
{
    Test::variantA();

    const auto &r = Standard::instance().r_var;

    REQUIRE(r.fold("D_1_1_R")  == 1);
    REQUIRE(r.fold("D_1_2_R")  == 2);
    REQUIRE(r.fold("D_1_3_R")  == 4);
    REQUIRE(r.fold("D_1_4_R")  == 8);
    REQUIRE(r.fold("D_1_5_R")  == 16);
    REQUIRE(r.fold("D_1_6_R")  == 32);
    REQUIRE(r.fold("D_1_7_R")  == 64);
    REQUIRE(r.fold("D_1_8_R")  == 128);
    REQUIRE(r.fold("D_1_9_R")  == 256);
    REQUIRE(r.fold("D_1_10_R") == 512);
    REQUIRE(r.fold("D_1_11_R") == 1019);
    REQUIRE(r.fold("D_1_12_R") == 2040);

    REQUIRE(r.fold("D_2_1_R")  == 1);
    REQUIRE(r.fold("D_2_2_R")  == 2);
    REQUIRE(r.fold("D_2_3_R")  == 4);
    REQUIRE(r.fold("D_2_4_R")  == 8);
    REQUIRE(r.fold("D_2_5_R")  == 16);
    REQUIRE(r.fold("D_2_6_R")  == 32);
    REQUIRE(r.fold("D_2_7_R")  == 64);
    REQUIRE(r.fold("D_2_8_R")  == 128);
    REQUIRE(r.fold("D_2_9_R")  == 256);
    REQUIRE(r.fold("D_2_10_R") == 512);
    REQUIRE(r.fold("D_2_11_R") == 1019);
    REQUIRE(r.fold("D_2_12_R") == 2040);

    REQUIRE(r.fold("D_3_1_R")  == 1);
    REQUIRE(r.fold("D_3_2_R")  == 2);
    REQUIRE(r.fold("D_3_3_R")  == 4);
    REQUIRE(r.fold("D_3_4_R")  == 8);
    REQUIRE(r.fold("D_3_5_R")  == 16);
    REQUIRE(r.fold("D_3_6_R")  == 32);
    REQUIRE(r.fold("D_3_7_R")  == 64);
    REQUIRE(r.fold("D_3_8_R")  == 128);
    REQUIRE(r.fold("D_3_9_R")  == 256);
    REQUIRE(r.fold("D_3_10_R") == 512);
    REQUIRE(r.fold("D_3_11_R") == 1019);
    REQUIRE(r.fold("D_3_12_R") == 2040);
}