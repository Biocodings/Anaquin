#ifndef SUBSAMPLER_HPP
#define SUBSAMPLER_HPP

#include <klib/khash.h>

namespace Anaquin
{
    class Subsampler
    {
        public:
            Subsampler(const std::string &hash, double prob) : _hash(hash), _prob(prob)
            {
                _seed = rand();
            }

            inline bool select() const
            {
                const uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(_hash.c_str()) ^ _seed);
                return ((double)(k&0xffffff) / 0x1000000 >= _prob);
            }

        private:
        
            const std::string _hash;
        
            // Random seed
            int _seed;

            // The probability of selection
            const double _prob;
    };
}

#endif