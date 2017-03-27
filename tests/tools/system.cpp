#include <fstream>
#include <catch.hpp>
#include "tools/system.hpp"

using namespace Anaquin;

TEST_CASE("Path2file_1")
{
    REQUIRE(path2file("/tmp/filegGPbYKhCDmQB/ABCD.gtf") == "ABCD.gtf");
    REQUIRE(path2file("/filegGPbYKhCDmQB/ThisIsTest.txt") == "ThisIsTest.txt");
}

TEST_CASE("ScriptTool_1")
{
    std::ifstream file("tests/data/Invalid.R");
    std::string str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    
    const auto r = System::trim(str);
    
    REQUIRE(r.size() == 780);
}
