#pragma once

#include <QWidget>
#include <QPainter>
#include <QColor>
#include "utils.h"

class Ploter;

class Axis : public QWidget
{
	Q_OBJECT
public:
	explicit Axis(Ploter & ploter, QColor color,
				  QWidget *parent = 0);

	void paintEvent(QPaintEvent * event);

	void drawOX(QPainter & p);
	void drawOY(QPainter & p);
signals:

public slots:

private:
	Ploter & m_ploter;
	QColor m_color;
	double m_min_point;

	static const int MIN_WIDTH;
	static const int MAJ_WIDTH;
};

#include "ploter.h"
