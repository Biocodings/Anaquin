#include <fstream>
#include <iostream>
#include <catch.hpp>

extern int parse_options(const std::string &command, std::string &output, std::string &error);

TEST_CASE("Meta_Print_Sequins")
{
    std::string output, error;
    const int status = parse_options("meta -l", output, error);

    REQUIRE(status == 0);
    REQUIRE(output.find("MG_30") != std::string::npos);
    REQUIRE(output.find("MG_29") != std::string::npos);
    REQUIRE(output.find("MG_23") != std::string::npos);
    REQUIRE(output.find("MG_47") != std::string::npos);
    REQUIRE(output.find("Metagenomics command detected") != std::string::npos);
}

TEST_CASE("Meta_Invalid_Filters")
{
    std::string output, error;

    std::ofstream o;
    o.open ("test.filter");
    o << "M11_G" << std::endl;
    o << "This is my sequin!!!" << std::endl;
    o.close();

    const int status = parse_options("meta -assembly contig.fa -r test.filter", output, error);
    
    REQUIRE(status == 1);
    REQUIRE(error.find("Unknown") != std::string::npos);
    REQUIRE(error.find("This is my sequin!!!") != std::string::npos);

    remove("test.filer");
}