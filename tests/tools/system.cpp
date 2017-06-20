#include <fstream>
#include <catch.hpp>
#include "tools/system.hpp"

using namespace Anaquin;

TEST_CASE("ScriptTool_1")
{
    std::ifstream file("tests/data/Invalid.R");
    std::string str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    
    const auto r = System::trim(str);
    
    REQUIRE(r.size() == 780);
}
