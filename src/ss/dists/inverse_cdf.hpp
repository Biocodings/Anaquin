#ifndef SS_INVERSE_CDF_HPP
#define SS_INVERSE_CDF_HPP

namespace SS
{
    template <typename P, typename I, typename T> struct InverseCDF
    {
        T operator()()
        {
            return i(p.next());
        }

        P p;
        I i;
    };

	/*
	 * Inverse cumulative normal distribution
	 *
	 * Given a value between zero and one as the cumulative probability
	 * of a gaussian normal distribution this class provides the value
	 * that it produces the specified CDF.
	 *
	 * For example, CDF(0) = 0.5 for N(0,1), so if we give 0.5 the result
	 * will be 0.
	 *
     * It use Acklam's approximation:
     * by Peter J. Acklam, University of Oslo, Statistics Division.
     * URL: http://home.online.no/~pjacklam/notes/invnorm/index.html
	 *
	 * This class is an extract from QuantLib.
	 */

	class InverseCumulativeNormal
	{
		public:
			InverseCumulativeNormal(double average = 0.0, double sigma = 1.0);
			
			double operator()(double x) const
			{
				return average_ + sigma_*standard_value(x);
			}

			/*
			 * Compared to operator(), this method avoids 2 floating point
			 * operations (we use average=0 and sigma=1 most of the
			 * time). The speed difference is noticeable.
	  		 */

			static double standard_value(double x)
			{
				double z;

				if (x < x_low_ || x_high_ < x)
				{
					z = tail_value(x);
				}
				else
				{
					z = x - 0.5;
					double r = z*z;
					z = (((((a1_*r + a2_)*r + a3_)*r + a4_)*r + a5_)*r + a6_)*z /
						(((((b1_*r + b2_)*r + b3_)*r + b4_)*r + b5_)*r + 1.0);
				}

				// The relative error of the approximation has absolute value less
				// than 1.15e-9.  One iteration of Halley's rational method (third
				// order) gives full machine precision.
				// #define REFINE_TO_FULL_MACHINE_PRECISION_USING_HALLEYS_METHOD

				return z;
			}

		private:
			static double tail_value(double x);
			//static const CumulativeNormalDistribution f_;

			double average_, sigma_;
			static const double a1_;
			static const double a2_;
			static const double a3_;
			static const double a4_;
			static const double a5_;
			static const double a6_;
			static const double b1_;
			static const double b2_;
			static const double b3_;
			static const double b4_;
			static const double b5_;
			static const double c1_;
			static const double c2_;
			static const double c3_;
			static const double c4_;
			static const double c5_;
			static const double c6_;
			static const double d1_;
			static const double d2_;
			static const double d3_;
			static const double d4_;
			static const double x_low_;
			static const double x_high_;
	};

	inline InverseCumulativeNormal::InverseCumulativeNormal(double average, double sigma)
		: average_(average), sigma_(sigma) {
		/*
		QL_REQUIRE(sigma_>0.0,
			"sigma must be greater than 0.0 ("
			<< sigma_ << " not allowed)");
			*/
	}
}

#endif