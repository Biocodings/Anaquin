//#include <catch.hpp>
//#include "rna/r_align.hpp"
//
//extern int parse_options(const std::string &, std::string &, std::string &);
//
//TEST_CASE("T_RNA_Cufflinks")
//{
//    std::string output, error;
//
//    const int s1 = parse_options("rna -align tests/data/rna/cufflinks.sam", output, error);
//    
//    REQUIRE(s1 == 0);
//    
//    const int s2 = parse_options("rna -model invalid_model -mixture invalid_mixture -align tests/data/rna/cufflinks.sam", output, error);
//
//    REQUIRE(s2 == 1);
//    REQUIRE(output.find("Model file")   != std::string::npos);
//    REQUIRE(output.find("Mixture file") != std::string::npos);
//}