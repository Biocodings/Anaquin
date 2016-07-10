#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <klib/khash.h>
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct Alignment;

    struct Sampler
    {
        struct SGReads
        {
            Reads syn = 0, gen = 0;
            
            inline Proportion dilut() const
            {
                return static_cast<Proportion>(syn) / (syn + gen);
            }
        };

        struct Stats
        {
            SGReads before, after;
        };
        
        struct User
        {
            virtual void syncReadSampled(const Alignment &) = 0;
        };

        static Stats subsample(const FileName &, Proportion, const AnalyzerOptions &, User *user = nullptr);
    };
    
    class Random
    {
        public:
        
            Random(double prob) : _prob(prob)
            {
                assert(prob >= 0.0);
                _seed = rand();
            }
        
            inline bool select(const std::string &hash) const
            {
                const uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(hash.c_str()) ^ _seed);
                return ((double)(k&0xffffff) / 0x1000000 >= _prob);
            }
        
        private:
        
            // Random seed
            int _seed;

            // The probability of selection
            Probability _prob;
    };
}

#endif