#include <fstream>
#include <streambuf>
#include <catch.hpp>
#include "tools/script.hpp"
#include <iostream>
using namespace Anaquin;

TEST_CASE("ScriptTool_1")
{
    std::ifstream file("tests/data/Invalid.R");
    std::string str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

    const auto r = Script::trim(str);

    REQUIRE(r.size() == 780);
}