#include "ploter.h"

Ploter::Ploter(double x_min, double x_max,
			   double y_min, double y_max,
               QColor axis_color,
               const char * title)
	: QWidget(nullptr)
	, m_x_min(x_min)
	, m_x_max(x_max)
	, m_y_min(y_min)
	, m_y_max(y_max)
	, m_axis(*this, axis_color, this)
{
    setWindowTitle(tr(title));
	resize(QSize(640, 480));
}

void Ploter::paintEvent(QPaintEvent *)
{
	QPainter p(this);
	for(auto plot: m_plots)
	{
		p.setPen(plot->get_color());
		for (double x = plot->x_min(),
					max = plot->x_max(),
					step = 1. / m_scale;
			 x < max;
			 x += step)
		{
			p.drawLine(screen_coord(x, plot->get_val(x)),
					   screen_coord(
						   x + step,
						   plot->get_val(x + step)
					   ));
		}
	}
}

void Ploter::resizeEvent(QResizeEvent *)
{
	m_scale = qMin(
		(geometry().width() / (m_x_max - m_x_min)),
		(geometry().height() / (m_y_max - m_y_min))
	);
	/*
	m_null_x = qMax(qFloor(m_scale * ( -m_x_min )),
					geometry().width() / 2);
	m_null_y = qMax(qFloor(m_scale * ( m_y_max )),
					geometry().height() / 2);
	*/
	m_null_x = qFloor(m_scale * ( -m_x_min ));
	m_null_y = qFloor(m_scale * ( m_y_max ));
	m_axis.resize(size());
}

void Ploter::drawPlot(Plot * plot)
{
	m_plots.append(plot);
	update();
}

QPoint Ploter::screen_coord(double x, double y)
{
	return QPoint(( (x * m_scale) + m_null_x),
				  (-(y * m_scale) + m_null_y));
}
