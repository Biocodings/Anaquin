#include <catch.hpp>
#include "tools/sample.hpp"

using namespace Anaquin;

TEST_CASE("Random_None")
{
    Random r(1.0);
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (r.select(std::to_string(100 * i))) { n++; }
    }
    
    REQUIRE(n == 0);
}

TEST_CASE("Random_Half")
{
    Random r(0.5);
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (r.select(std::to_string(100 * i))) { n++; }
    }

    REQUIRE(n > 10);
    REQUIRE(n < 90);
}

TEST_CASE("SamplingTool_All")
{
    Random r(0.0);
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (r.select(std::to_string(100 * i))) { n++; }
    }
    
    REQUIRE(n == 100);
}