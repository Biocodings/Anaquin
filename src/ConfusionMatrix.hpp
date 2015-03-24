#ifndef AS_CONFUSION_MATRIX_HPP
#define AS_CONFUSION_MATRIX_HPP

#include <math.h>
#include "Types.hpp"

struct ConfusionMatrix
{
	Percentage fp = 0;
	Percentage tp = 0;
	Percentage fn = 0;
	Percentage tn = 0;

	inline Percentage sp() const
	{
		return (tp + fn) ? tp / (tp + fn) : NAN;
	}

	inline Percentage sn() const
	{
		return (fp + tn) ? tn / (fp + tn) : NAN;
	}
};

#endif