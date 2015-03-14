#include "gtest/gtest.h"
#include "../AlignerCompare.hpp"

using namespace std;

TEST(CufflinkExample, AlignerStatsTest)
{
	/*
	 * The sample file was taken from Cufflink's source code. It's obviously independent to the standards.
	 */

	const auto stats = AlignerCompare::analyze("C:\\Sources\\QA\\Tests\\Data\\CufflinksTest.sam");

	/*
	 * There shouldn't be any match to the in-silico chromosome. Sensivitiy is zero because our experiment
	 * has failed to detect anything from the standards. Specificity should be one because none of the reads
	 * that fails to map to the chromosome comes from it.
	 */



}