#ifndef SS_LEVENE_TEST_HPP
#define SS_LEVENE_TEST_HPP

#include "QStatistics/Distribution/FDist.hpp"
#include "QStatistics/tests/StatisticalTest.hpp"

namespace SS
{
	/*
	 * Levene's test is a statistic test used to assess the equality of variances
	 * among two or more groups. Some common statistical algorithms assume that variances
	 * of the populations from which different samples are drawn are equal. Levene's test
	 * assesses this assumption. It tests the null hypothesis that the population variances
	 * are equal, also known as homoscedasticity.
	 *
	 * The null hypothesis is that variances of all groups are equal. The alternative hypothesis
	 * is that at least one of the variances is significantly different. Levene's test is always
	 * one-tailed, and uses an F statistic.
	 */

	/*
	 * A matlab implementation can be found at
	 *      http://au.mathworks.com/matlabcentral/fileexchange/3375-levenetest/content/Levenetest.m
	 */

	template <typename Iter, typename T> Results leveneTest(const Iter &iter)
	{
		static_assert(std::is_same<NumericalVariable, Iter::value_type>::value, "NumericalVariable is expected");

		// Total number of cases in all groups
		SampleSize n = 0;

		// Number of different groups
		const auto k = iter.size();

		// A matrix of transformed data (columns needn't have the same size)
		std::vector<std::vector<T>> Z(k);

		// Mean for each group
		std::vector<T> uk(k);

		// Running sum, needed for the overall mean
		T r_sum = 0;

		/*
		 * Loop over each group and compute the statistics
		 */

		std::size_t j = 0;

		std::for_each(iter.begin(), iter.end(), [&](const NumericalVariable &v)
		{
			const auto u = mean(v.values());

			// Running sum for a group
			T sum = 0;

			for (std::size_t i = 0; i < v.values().size(); i++)
			{
				n++;

				// Transform a raw value to a distance from the mean or medium
				const auto z = fabs(v.values()[i] - u);

				// The individual values will be needed for computing a chi-square
				Z[j].push_back(z);

				// Update the running sum for this group
				sum += z;
			};

			// Update the overall running sum
			r_sum += sum;

			// Compute the mean for this group
			uk[j++] = sum / v.values().size();
		});

		// The overall mean
		const auto u = r_sum / n;

		// Degree of freedom for the numerator
		const DF df_n = k - 1;

		// Degree of freedom for the denominator
		const DF df_d = n - k;

		// Chi-square for the numerator with df_n degree of freedom
		const auto num = sum<T>(0, k, [&](std::size_t i)
		{
			return Z[i].size() * std::pow((uk[i] - u), 2);
		});

		// Chi-square for the denominator with df_d degree of freedom
		const auto den = sum<T>(0, k, [&](std::size_t j)
		{
			return sum<T>(0, Z[j].size(), [&](std::size_t i)
			{
				return (Z[j][i] - uk[j]) * (Z[j][i] - uk[j]);
			});
		});

		TestResults<T> r;

		// The statistic is a F-statistic of two independent chi-squares
		r.x = (num * df_d) / (den * df_n);

		r.pv = FDist_T<T>(df_n, df_d).cdf(r.x);

		return r;
	}

	/*
	 * Perform a levene test between two variables
	 */

	template<typename T> TestResults<T> leveneTest(const NumericalVariable &v1, const NumericalVariable &v2)
	{
		const std::vector<NumericalVariable> vs = { v1, v2 };
		return leveneTest<std::vector<NumericalVariable>, T>(vs);
	}
}

#endif