#ifndef SS_MERSENNE_TWISTER_HPP
#define SS_MERSENNE_TWISTER_HPP

#include <vector>

namespace SS
{
	/*
	 * Uniform random number generator.
	 * Mersenne Twister random number generator of period 2**19937-1.
	 *
	 *		For more details see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
	 *
	 * This file is an extract from QuantLib.
	 */

    class MersenneTwisterUniformRng
	{
      public:
        //typedef Sample<double> sample_type;

        /*! if the given seed is 0, a random seed will be chosen
            based on clock() */
        explicit MersenneTwisterUniformRng(unsigned long seed = 0);
        explicit MersenneTwisterUniformRng(
                                     const std::vector<unsigned long>& seeds);
        /*! returns a sample with weight 1.0 containing a random number
            in the (0.0, 1.0) interval  */
		double next() const { return nextReal(); }
        //! return a random number in the (0.0, 1.0)-interval
		double nextReal() const {
			return (double(nextInt32()) + 0.5) / 4294967296.0;
        }
        //! return a random int in the [0,0xffffffff]-interval
        unsigned long nextInt32() const  {
            if (mti==N)
                twist(); /* generate N words at a time */

            unsigned long y = mt[mti++];

            /* Tempering */
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680UL;
            y ^= (y << 15) & 0xefc60000UL;
            y ^= (y >> 18);
            return y;
        }

      private:
		static const unsigned N = 624; // state size
		static const unsigned M = 397; // shift size
        void seedInitialization(unsigned long seed);
        void twist() const;
        mutable unsigned long mt[N];
        mutable unsigned mti;
        static const unsigned long MATRIX_A, UPPER_MASK, LOWER_MASK;
    };
}

#endif