#include "function.h"

double Sum::value(std::map<char, double> var_vals)
{
	return (m_left->value(var_vals) +
			m_right->value(var_vals));
}

Function * Sum::derivative(char var)
{
	return new Sum(m_left->derivative(var),
				   m_right->derivative(var));
}

