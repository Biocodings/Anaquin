#include <catch.hpp>

extern int parse_options(const std::string &command, std::string &output);

TEST_CASE("Meta_Print_Sequins")
{
    std::string output;
    const int status = parse_options("meta -l", output);

    REQUIRE(status == 0);
    REQUIRE(output.find("MG_30") != std::string::npos);
    REQUIRE(output.find("MG_29") != std::string::npos);
    REQUIRE(output.find("MG_23") != std::string::npos);
    REQUIRE(output.find("MG_47") != std::string::npos);
    REQUIRE(output.find("Metagenomics command detected") != std::string::npos);
}