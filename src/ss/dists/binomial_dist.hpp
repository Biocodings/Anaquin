#ifndef QQ_BINOMIAL_DIST_HPP
#define QQ_BINOMIAL_DIST_HPP

#include "QStatistics/Statistics.hpp"
#include "QStatistics/Distribution/Distribution.hpp"

namespace QQ
{
	/*
	 * The binomial distribution is used to model the number of successes in a sample of size n drawn
	 * with replacement from a population of size N.
	 */

	template<typename T> struct BinomialDist_T : public Distriubtion_<T>
	{
		BinomialDist_T(Count n, Probability p) : n(n), p(p) {}

		T mean() const;
		T variance() const;

		// The size of the sample
		Count n;

		// The probability of success for each sample
		Probability p;
	};

	template<typename T> T BinomialDist_T<T>::mean() const
	{
		return T(n * p);
	}

	template<typename T> T BinomialDist_T<T>::variance() const
	{
		return T(n * p * (1 - p));
	}

	typedef BinomialDist_T<Real> BinomialDist;
}

#endif