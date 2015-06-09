#include <fstream>
#include <iostream>
#include <catch.hpp>

extern int parse_options(const std::string &command, std::string &output, std::string &error);

TEST_CASE("Meta_Single_Filter")
{
    std::string output, error;
    
    std::ofstream o;
    o.open ("test.filter");
    o << "M11_G" << std::endl;
    o.close();
    
    const int status = parse_options("meta -assembly tests/data/meta/e1/contigs_A.fa -f test.filter -psl tests/data/meta/e1/align_A.psl", output, error);

    REQUIRE(status == 0);
    REQUIRE(output.find("Metagenomics") != std::string::npos);
    
    remove("test.filter");
}

TEST_CASE("Meta_Invalid_Filters")
{
    std::string output, error;
    
    std::ofstream o;
    o.open ("test.filter");
    o << "M11_G" << std::endl;
    o << "This is my sequin!!!" << std::endl;
    o.close();
    
    const int status = parse_options("meta -assembly contig.fa -f test.filter", output, error);
    
    REQUIRE(status == 1);
    REQUIRE(error.find("Unknown") != std::string::npos);
    REQUIRE(error.find("This is my sequin!!!") != std::string::npos);
    
    remove("test.filter");
}

TEST_CASE("Meta_Assembly_E1")
{
    std::string output, error;

    const int status_1 = parse_options("meta ./anaquin meta -assembly tests/data/meta/e1/contigs_A.fa -psl tests/data/meta/e1/align_A.psl", output, error);

    REQUIRE(status_1 == 0);

    const int status_2 = parse_options("meta ./anaquin meta -diffs tests/data/meta/e1/contigs_A.fa tests/data/meta/e1/contigs_B.fa tests/data/meta/e1/align_A.psl tests/data/meta/e1/align_B.psl", output, error);

    REQUIRE(status_2 == 0);
}

///*
// * Analyze an assembly from Velvet. Given the same default mixture and annotation model.
// */
//
//TEST_CASE("Meta_Assembly_Default_Mixture_Model")
//{
//    std::string output, error;
//    
//    const int status = parse_options("meta -mixture data/meta/META.mix.csv -model data/meta/META.ref.bed -assembly tests/data/meta/contigs_3.fa", output, error);
//
//    REQUIRE(status == 0);
//    REQUIRE(output.find("Mixture file") != std::string::npos);
//    REQUIRE(output.find("Model file") != std::string::npos);
//}

///*
// * Analyze an assembly from Velvet
// */
//
//TEST_CASE("Meta_Assembly")
//{
//    std::string output, error;
//    
//    const int status = parse_options("meta -assembly tests/data/meta/contigs_3.fa", output, error);
//    
//    REQUIRE(status == 0);
//}

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

TEST_CASE("Meta_Blast")
{
//    std::string output, error;
    
//    const int status = parse_options("meta -blast tests/data/meta/align.psl", output, error);

//    REQUIRE(status == 0);
//    REQUIRE(output.find("M5_G") != std::string::npos);
}