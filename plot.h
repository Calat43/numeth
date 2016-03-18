#pragma once

#include <QColor>

class Plot
{
public:
	Plot(double x_min, double x_max,
		 double y_min, double y_max,
		 QColor color);
	Plot(const Plot & src);

	virtual double get_val(double x) = 0;

	QColor & get_color()
	{
		return m_color;
	}
	double x_min() const
	{
		return m_x_min;
	}
	double x_max() const
	{
		return m_x_max;
	}
	double y_min() const
	{
		return m_y_min;
	}
	double y_max() const
	{
		return m_y_max;
	}

private:
	double m_x_min;
	double m_x_max;
	double m_y_min;
	double m_y_max;
	QColor m_color;
};
