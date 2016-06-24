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

            inline Counts &tp() const { return _tp; }
            inline Counts &fp() const { return _fp; }
            inline Counts &tn() const { return _tn; }
            inline Counts &fn() const { return _fn; }
            inline Counts &nr() const { return _nr; }
            inline Counts &nq() const { return _nq; }

            // Calculate for the sensitivity
            inline Proportion sn() const
            {
                assert(_nr && _nr >= _tp);

                // Adjust for fn... Refer to wikipedia for details
                _fn = _nr - _tp;

                return (_tp + _fn) ? static_cast<Proportion>(_tp) / (_tp + _fn) : NAN;
            }

            // Calculate for the specificity
            inline Proportion sp() const
            {
                assert(_nr && _nr >= _tp);
                return (_tn + _fp) ? static_cast<Proportion>(_tn) / (_tn + _fp) : NAN;
            }

            // Calculate for the precision
            inline Proportion pc() const
            {
                assert(_nr && _nr >= _tp);
                return (_tp + _fp) ? static_cast<Proportion>(_tp) / (_tp + _fp) : NAN;
            }

            // Calculate for the FDR
            inline Proportion fdr() const
            {
                return 1.0 - pc();
            }

        private:
        
            mutable Counts _fp = 0;
            mutable Counts _tp = 0;
            mutable Counts _fn = 0;
            mutable Counts _tn = 0;
            mutable Counts _nq = 0;
            mutable Counts _nr = 0;
    };
}

#endif