#include <fstream>
#include <iostream>
#include <catch.hpp>

extern int parse_options(const std::string &command, std::string &output, std::string &error);

// Generate a file for filtering, one line per sequin
static void generateFilter(const std::string &file, const std::vector<std::string> &sequins)
{
    std::ofstream o;
    o.open(file);
    
    for (const auto &sequin : sequins)
    {
        o << sequin << std::endl;
    }
    
    o.close();
}

TEST_CASE("Meta_Assembly")
{
    std::string output, error;
    
    const int s1 = parse_options("meta -assembly tests/data/meta/A/contigs.fa -psl tests/data/meta/A/align.psl", output, error);
    REQUIRE(s1 == 0);

    const int s2 = parse_options("meta -diffs tests/data/meta/A/contigs.fa tests/data/meta/B/contigs.fa -psl tests/data/meta/A/align.psl tests/data/meta/B/align.psl", output, error);
    REQUIRE(s2 == 0);
}

TEST_CASE("Meta_Diffs_Empty_Align")
{
    std::string output, error;

    // Typical command for differential with alignments
    const int s = parse_options("meta -diffs tests/data/meta/A/contigs.fa tests/data/meta/B/contigs.fa -psl tests/data/meta/empty.psl tests/data/meta/empty.psl", output, error);

    REQUIRE(s == 1);
    REQUIRE(error.find("Empty file:") != std::string::npos);
}

TEST_CASE("Meta_Sequins")
{
    std::string output, error;
    const int status = parse_options("meta -l", output, error);
    
    REQUIRE(status == 0);
    REQUIRE(output.find("MG_30") != std::string::npos);
    REQUIRE(output.find("MG_29") != std::string::npos);
    REQUIRE(output.find("MG_23") != std::string::npos);
    REQUIRE(output.find("MG_47") != std::string::npos);
    REQUIRE(output.find("Metagenomics") != std::string::npos);
}

TEST_CASE("Meta_Single_Filter")
{
    std::string output, error;
    generateFilter("test.filter", std::vector<std::string> { "M11_G" });
    
    const int status = parse_options("meta -assembly tests/data/meta/A/contigs.fa -f test.filter -psl tests/data/meta/A/align.psl", output, error);

    REQUIRE(status == 0);
    REQUIRE(output.find("Metagenomics") != std::string::npos);
    
    remove("test.filter");
}

TEST_CASE("Meta_Invalid_Filters")
{
    std::string output, error;
    generateFilter("test.filter", std::vector<std::string> { "M11_G", "This is my sequin!!!" });
    
    const int status = parse_options("meta -assembly contig.fa -f test.filter", output, error);
    
    REQUIRE(status == 1);
    REQUIRE(error.find("Unknown") != std::string::npos);
    REQUIRE(error.find("This is my sequin!!!") != std::string::npos);
    
    remove("test.filter");
}