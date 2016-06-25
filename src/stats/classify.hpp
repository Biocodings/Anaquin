#ifndef CLASSIFY_HPP
#define CLASSIFY_HPP

#include <cmath>
#include <assert.h>
#include "data/types.hpp"

namespace Anaquin
{
    class Confusion
    {
        public:

            void operator+=(const Confusion &m)
            {
                _tp += m.tp();
                _fp += m.fp();
                _fn += m.fn();
                _nr += m.nr();
                _nq += m.nq();
            }

            inline Counts &tp() const { return _tp; }
            inline Counts &fp() const { return _fp; }
            inline Counts &fn() const { return _fn; }
            inline Counts &nr() const { return _nr; }
            inline Counts &nq() const { return _nq; }

            inline Proportion sn() const
            {
                return (_tp + _fn) ? static_cast<Proportion>(_tp) / (_tp + _fn) : NAN;
            }

            inline Proportion pc() const
            {
                return (_tp + _fp) ? static_cast<Proportion>(_tp) / (_tp + _fp) : NAN;
            }

            inline Proportion fdr() const
            {
                return 1.0 - pc();
            }

        private:
        
            mutable Counts _fp = 0;
            mutable Counts _tp = 0;
            mutable Counts _fn = 0;
            mutable Counts _nq = 0;
            mutable Counts _nr = 0;
    };
}

#endif