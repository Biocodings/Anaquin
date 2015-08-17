#include <catch.hpp>
#include "data/standard.hpp"

using namespace Anaquin;

TEST_CASE("Standard_ChrT")
{
	REQUIRE(Standard::instance().id == "chrT");
}