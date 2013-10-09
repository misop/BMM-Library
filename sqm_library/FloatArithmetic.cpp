#include "stdafx.h"
#include "FloatArithmetic.h"
#include <math.h>

bool equal(float x, float y) {
	return fabs(x - y) < FLOAT_ZERO;
}