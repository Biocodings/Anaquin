#include <fstream>
#include <streambuf>
#include <catch.hpp>
#include "tools/script.hpp"

using namespace Anaquin;

TEST_CASE("ScriptTool_1")
{
    std::ifstream file("tests/data/Invalid.R");
    std::string str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    const auto r = ScriptTool::clean(str);

    REQUIRE(r.size() == 754);
}