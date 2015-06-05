#include <catch.hpp>

extern int parse_options(const std::string &command, std::string &output);

TEST_CASE("Meta_Print_Sequins")
{
    std::string output;
    const int status = parse_options("meta -l", output);

    std::cout << output << std::endl;
    
    REQUIRE(status == 0);
    
}

