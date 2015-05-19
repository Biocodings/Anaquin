#ifndef SS_CODING_HPP
#define SS_CODING_HPP

#include <set>
#include <ss/matrix.hpp>

namespace SS
{
	struct DummyCoding
	{
		template<typename I1, typename I2> static Matrix encode(const I1 &y, const I2 &f);
	};

	template <typename I1, typename I2> Matrix DummyCoding::encode(const I1 &y, const I2 &f)
	{
		// Number of regression coefficients
		const auto p = std::set<int>(f.begin(), f.end()).size() - 1;

		// The regression coefficient matrix (a column for the constant coefficient)
		Matrix X(y.size(), p + 1);

        X.setConstant(0);

        // Constant coefficient
        X.col(0).setOnes();

		for (auto i = 0; i < y.size(); i++)
		{
			auto code = (int)f[i] > 1 ? 1 << ((int)f[i] - 2) : 0;

			for (auto j = p; j && code; j--, code >>= 1)
			{
				if (code & 1)
				{
					X(i, j) = 1;
					break;
				}
			}
		}

		return X;
	}
}

#endif