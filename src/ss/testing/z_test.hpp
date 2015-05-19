#ifndef SS_Z_TEST_HPP
#define SS_Z_TEST_HPP

#include <sstats/dists/normal.hpp>
#include "QStatistics/Testing/StatisticalTest.hpp"

namespace QQ
{
	template<typename T> TestResults<T> oneSampleZTest(SampleSize n, T u, T s, StatisticalTestType type, T u0 = 0)
	{
		return statTest<T, TestResults<T>>((u - u0) / (s / std::sqrt(n)), StandardNormal_T<T>(), type);
	}
}

#endif