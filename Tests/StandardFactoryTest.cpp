#include "gtest/gtest.h"
#include "StandardFactory.hpp"

using namespace std;

TEST(TestChromoName, StandardFactoryTest)
{
	const auto name = StandardFactory::chromoName();
	ASSERT_EQ(name, "chrT");
}