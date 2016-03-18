#include "utils.h"

bool utils::is_int(double val, double eps)
{
	return (floor(val + 0.5) - eps <= val &&
			val <= floor(val + 0.5) + eps);
}

