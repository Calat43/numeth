#pragma once

#include <set>
#include <map>

class Function
{
public:
	//Function();
	//Function(const Function & src);
	virtual ~Function() {}

	virtual double
		value(std::map<char, double> var_vals) = 0;
	virtual Function * derivative(char var) = 0;
};

class Sum : public Function
{
public:
	Sum(Function * left, Function * right)
		: m_left(left)
		, m_right(right)
	{}
	Sum(const Sum & src)
		: m_left(src.m_left)
		, m_right(src.m_right)
	{}
	~Sum()
	{
		delete m_left;
		delete m_right;
	}
	virtual double
		value(std::map<char, double> var_vals);
	virtual Function * derivative(char var);
private:
	Function * m_left;
	Function * m_right;
	std::set<char> m_vars;
};


