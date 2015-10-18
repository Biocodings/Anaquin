#include <catch.hpp>
#include "tools/sampling.hpp"

using namespace Anaquin;

TEST_CASE("SamplingTool_None")
{
    SamplingTool sampler(1.0);
    
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (sampler.select(std::to_string(100 * i))) { n++; }
    }
    
    REQUIRE(n == 0);
}

TEST_CASE("SamplingTool_Half")
{
    SamplingTool sampler(0.5);
    
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (sampler.select(std::to_string(100 * i))) { n++; }
    }

    REQUIRE(n > 10);
    REQUIRE(n < 90);
}

TEST_CASE("SamplingTool_All")
{
    SamplingTool sampler(0.0);
    
    std::size_t n = 0;
    
    for (auto i = 0; i < 100; i++)
    {
        if (sampler.select(std::to_string(100 * i))) { n++; }
    }
    
    REQUIRE(n == 100);
}