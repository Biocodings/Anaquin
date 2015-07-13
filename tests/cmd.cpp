#include <catch.hpp>

extern int parse_options(const std::string &, std::string &, std::string &);

static std::string _output, _error;

TEST_CASE("Ladder_Sequis")
{
    const int s = parse_options("-c ladder -p seqs -m data/ladder/CON.v3.mix.csv", _output, _error);
    
    REQUIRE(s == 0);
}

TEST_CASE("Missing_Mixture_And_Reference")
{
    const int s = parse_options("-c ladder -p correct tests/data/ladder/aligned_A.sam", _output, _error);

    REQUIRE(s == 1);
    REQUIRE(_output.find("Ladder Analysis") != std::string::npos);
    REQUIRE(_error.find("Mixture file is missing. Please specify it with -m.") != std::string::npos);
}

TEST_CASE("RNA_Missing_Reference")
{
    const int s = parse_options("-c rna -p align -m data/rna/RNA.v4.1.mix tests/data/rna/cufflinks.sam", _output, _error);

    REQUIRE(s == 1);
    REQUIRE(_output.find("RNA Analysis") != std::string::npos);
    REQUIRE(_error.find("Reference file is missing. Please specify it with -r.") != std::string::npos);
}

TEST_CASE("Ladder_Missing_Reference")
{
    /*
     * Ladder analysis doesn't require a reference
     */
    
    const int s = parse_options("-c ladder -p correct -m data/ladder/CON.v3.mix.csv tests/data/ladder/aligned_A.sam", _output, _error);

    REQUIRE(s == 1);
    REQUIRE(_output.find("Ladder Analysis") != std::string::npos);
}


