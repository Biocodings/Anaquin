#include <catch.hpp>
//
extern int parse_options(const std::string &, std::string &, std::string &);

static std::string _output, _error;

TEST_CASE("Test_Version")
{
    const int s = parse_options("-v", _output, _error);
    
    REQUIRE(s == 0);
    REQUIRE(_output == "Anaquin v1.1.01\n");
}
//
///*
//TEST_CASE("Missing_Mixture_And_Reference")
//{
//    const int s = parse_options("-c ladder -p abund tests/data/ladder/aligned_A.sam", _output, _error);
//
//    REQUIRE(s == 1);
//    REQUIRE(_output.find("Ladder Analysis") != std::string::npos);
//    REQUIRE(_error.find("Mixture file is missing. Please specify it with -m.") != std::string::npos);
//}
//*/
//
///*
//TEST_CASE("RNA_Missing_Reference")
//{
//    const int s = parse_options("-c rna -p align -m data/rna/RNA.v4.1.mix tests/data/rna/cufflinks.sam", _output, _error);
//
//    REQUIRE(s == 1);
//    REQUIRE(_output.find("RNA Analysis") != std::string::npos);
//    REQUIRE(_error.find("Reference file is missing. Please specify it with -r.") != std::string::npos);
//}
//*/
