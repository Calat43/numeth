#include "axis.h"

const int Axis::MIN_WIDTH = 2;
const int Axis::MAJ_WIDTH = 5;

Axis::Axis(Ploter & ploter, QColor color, QWidget *parent)
    : QWidget(parent)
    , m_ploter(ploter)
    , m_color(color)
    , m_min_point(0.04)
{
}

void Axis::drawOX(QPainter & p)
{
    int ox_top = m_ploter.null_y();
    if (ox_top < 0)
    {
        ox_top = 0;
    }
    else if (ox_top > geometry().height())
    {
        ox_top = geometry().height();
    }
    p.drawLine(
        QPoint(geometry().x(),     ox_top),
        QPoint(geometry().width(), ox_top)
    );

    for(double x = floor(m_ploter.m_x_min / m_min_point) * m_min_point;
        x <= m_ploter.m_x_max + m_min_point / 2;
        x += m_min_point)
    {
        int width;
        int x_screen = m_ploter.screen_coord(x, 0).x();

        if (utils::is_int(x, m_min_point / 2))
        {
            if (round(x) == 0)
            {
                width = 0;
            }
            else
            {
                width = MAJ_WIDTH;
                QString label = QString::number(round(x));
                p.drawText(QPoint(x_screen - width,
                                  ox_top + width*4),
                           label);
            }
        }
        else
        {
            width = MIN_WIDTH;
        }
        p.drawLine(
            QPoint(x_screen, ox_top - width),
            QPoint(x_screen, ox_top + width)
        );
    }
}


void Axis::drawOY(QPainter & p)
{
    int oy_left = m_ploter.null_x();
    if (oy_left < 0)
    {
        oy_left = 0;
    }
    else if (oy_left > geometry().width())
    {
        oy_left = geometry().width();
    }
    p.drawLine(
        QPoint(oy_left, geometry().y()),
        QPoint(oy_left, geometry().height())
    );

    for(double y = floor(m_ploter.m_y_min / m_min_point) * m_min_point;
        y <= m_ploter.m_y_max + m_min_point / 2;
        y += m_min_point)
    {
        int width;
        int y_screen = m_ploter.screen_coord(0, y).y();

        if (utils::is_int(y, m_min_point / 2))
        {
            if (round(y) == 0)
            {
                width = 0;
            }
            else
            {
                width = MAJ_WIDTH;
                QString label = QString::number(round(y));
                p.drawText(QPoint(oy_left - width*4,
                                  y_screen + 4),
                           label);
            }
        }
        else
        {
            width = MIN_WIDTH;
        }
        p.drawLine(
            QPoint(oy_left - width, y_screen),
            QPoint(oy_left + width, y_screen)
        );
    }
}

void Axis::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.setPen(m_color);

    drawOX(p);
    drawOY(p);
}
