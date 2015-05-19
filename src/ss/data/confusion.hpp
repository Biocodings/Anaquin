#ifndef SS_CONFUSION_HPP
#define SS_CONFUSION_HPP

#include <math.h>
#include <ss/types.hpp>

namespace SS
{
    class Confusion
    {
        public:

            inline Counts &tp() const { return _tp; }
            inline Counts &fp() const { return _fp; }
            inline Counts &tn() const { return _tn; }
            inline Counts &fn() const { return _fn; }

            inline Counts n() const
            {
                return _fp + _tp + _fn + _tn;
            }

            // Sensitivity, metrics for positive classification
            inline Percentage sn() const
            {
                return _tp && (_tp + _fn) ? static_cast<Percentage>(_tp) / (_tp + _fn) : NAN;
            }

            // Specificity, metrics for negative classification
            inline Percentage sp() const
            {
                return _tn && (_tn + _fp) ? static_cast<Percentage>(_tn) / (_tn + _fp) : NAN;
            }

        protected:
            mutable Counts _fp = 0;
            mutable Counts _tp = 0;
            mutable Counts _fn = 0;
            mutable Counts _tn = 0;
    };
}

#endif