#include <catch.hpp>
#include "tools/system.hpp"

using namespace Anaquin;

TEST_CASE("Path2file_1")
{
    REQUIRE(path2file("/tmp/filegGPbYKhCDmQB/ABCD.gtf") == "ABCD.gtf");
    REQUIRE(path2file("/filegGPbYKhCDmQB/ThisIsTest.txt") == "ThisIsTest.txt");
}