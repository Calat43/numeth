#pragma once

#include <QtGlobal>
#include <QtMath>
#include <QVector>
#include <QWidget>
#include <QPainter>
#include <QPoint>
#include "axis.h"
#include "plot.h"

//class Axis;
//class Plot;

class Ploter : public QWidget
{
	Q_OBJECT

	friend class Axis;
public:
	Ploter(double x_min, double x_max,
		   double y_min, double y_max,
		   QColor axis_color);

	void paintEvent(QPaintEvent *);

	void resizeEvent(QResizeEvent *);

	int null_x() const
	{
		return m_null_x;
	}
	int null_y() const
	{
		return m_null_y;
	}

	void drawPlot(Plot * plot);
private:
	QPoint screen_coord(double x, double y);
private:
	double m_x_min;
	double m_x_max;
	double m_y_min;
	double m_y_max;
	QVector<Plot *> m_plots;
	Axis m_axis;
	int m_null_x;
	int m_null_y;

	int m_scale;
};
