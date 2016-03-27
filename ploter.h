#pragma once

#include <QtGlobal>
#include <QtMath>
#include <QVector>
#include <QWidget>
#include <QPainter>
#include <QPoint>
#include <QTimer>
#include "axis.h"
#include "plot.h"
#include <vector>

//class Axis;
//class Plot;

class Ploter : public QWidget
{
    Q_OBJECT

    friend class Axis;
    friend class AnimPloter;
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

class AnimPloter : public Ploter
{
    Q_OBJECT
public:
    AnimPloter(double x_min, double x_max,
               double y_min, double y_max,
               QColor axis_color,
               unsigned int delay)
        : Ploter(x_min, x_max, y_min, y_max, axis_color)
        , m_cur_plot(0)
    {
        QTimer *timer = new QTimer(this);
        connect(timer, SIGNAL(timeout()), this, SLOT(update()));
        timer->start(delay);
    }
    void paintEvent(QPaintEvent *)
    {
        QPainter p(this);

        if (m_plots.empty()) return;

        bool has_new_frames = false;
        for (auto plot_anim : m_plots)
        {
            Plot * plot;
            if (m_cur_plot < plot_anim.size())
            {
                plot = plot_anim[m_cur_plot];
                has_new_frames = true;
            }
            else
            {
                plot = plot_anim[plot_anim.size()-1];
            }
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
        if (has_new_frames)
        {
            ++m_cur_plot;
        }
    }

    void drawAnimPlot(std::vector<Plot *> & plot_anim)
    {
        m_plots.push_back(plot_anim);
        update();
        m_cur_plot = 0;
    }
/*
public slots:
    void update()
    {
        ++m_cur_plot;
    }*/
private:
    std::vector<std::vector<Plot *>> m_plots;
    unsigned int m_cur_plot;
};
