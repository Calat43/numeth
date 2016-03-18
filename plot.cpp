#include "plot.h"

Plot::Plot(double x_min, double x_max,
		   double y_min, double y_max,
		   QColor color)
	: m_x_min(x_min)
	, m_x_max(x_max)
	, m_y_min(y_min)
	, m_y_max(y_max)
	, m_color(color)
{}

Plot::Plot(const Plot & src)
	: m_x_min(src.m_x_min)
	, m_x_max(src.m_x_max)
	, m_y_min(src.m_y_min)
	, m_y_max(src.m_y_max)
	, m_color(src.m_color)
{}

