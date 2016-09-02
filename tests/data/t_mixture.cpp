#include <catch.hpp>
#include "data/standard.hpp"

using namespace Anaquin;

TEST_CASE("Mixture_MRN040")
{
    REQUIRE_THROWS_AS(Standard::instance().addTMix(Reader("data/RnaQuin/MRN040_v001.csv")), std::runtime_error);
}

TEST_CASE("Mixture_MRN041")
{
    REQUIRE_THROWS_AS(Standard::instance().addTMix(Reader("data/RnaQuin/MRN041_v001.csv")), std::runtime_error);
}

TEST_CASE("Mixture_MRN042")
{
    REQUIRE_THROWS_AS(Standard::instance().addTMix(Reader("data/RnaQuin/MRN042_v001.csv")), std::runtime_error);
}