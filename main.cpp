#include "ploter.h"
#include <QApplication>
#include <QtMath>
#include "numeth2.h"
#include "numeth3.h"
#include <fstream>
/*
class Sin_plot : public Plot
{
public:
    Sin_plot(double x_min, double x_max,
             double y_min, double y_max,
             QColor color)
        : Plot(x_min, x_max, y_min, y_max, color) {}
    double get_val(double x)
    {
        return sin(x);
    }
};

class Array_plot : public Plot
{
public:
    Array_plot(double x_min, double x_max,
               double y_min, double y_max,
               QColor color,
               std::vector<double> data)
        : Plot(x_min, x_max, y_min, y_max, color)
        , m_data(data)
    {}

    double get_val(double x)
    {
        if (x < x_min() || x > x_max()) return 0;

        int id = (x - x_min()) /
                 ((x_max() - x_min()) / m_data.size());

        return m_data[id];
    }
private:
    std::vector<double> m_data;
};

class Array_plot2 : public Plot
{
public:
    Array_plot2(double x_min, double x_max,
               double y_min, double y_max,
               QColor color,
               std::vector<double> data)
        : Plot(x_min, x_max, y_min, y_max, color)
        , m_data(data)
    {}

    double get_val(double x)
    {
        if (x < x_min() || x > x_max()) return 0;

        int id = (x - x_min()) /
                 ((x_max() - x_min()) / m_data.size());

        //double h =
        return m_data[id];
    }
private:
    std::vector<double> m_data;
};*/

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
    std::ofstream fout("output.txt", std::ofstream::out);
    std::ofstream fout_appr("output_appr.txt", std::ofstream::out);
   // AnimPloter trian_ploter_win(-0.1, 2, -2, 2, QColor(0, 0, 0), 100, "triangle");
    //trian_ploter_win.show();

    AnimPloter ro_ploter_win(-4.1, 6.1, 0, 1, QColor(0, 0, 0), 100, "dencity");
    AnimPloter u_ploter_win(-4.1, 6.1, 0, 1, QColor(0, 0, 0), 100, "velocity");
    AnimPloter p_ploter_win(-4.1, 6.1, 0, 1, QColor(0, 0, 0), 100, "pressure");
    ro_ploter_win.show();
    u_ploter_win.show();
    p_ploter_win.show();
    //task1(ploter_win, fout);
    //expl(ploter_win, fout);
    //triangle(trian_ploter_win, fout_appr);
    maccormack(ro_ploter_win, u_ploter_win, p_ploter_win, fout);


    /*
    Ploter ploter_win(-2, 2, -3, 3, QColor(0,0,0));
    ploter_win.show();
    std::vector<double> plot_data;
    for (int i = 0; i < 20; i++)
    {
        plot_data.push_back(sin(-2 + (4. / 20) * i));
    }
    Plot * plot = new Array_plot(-2, 2, -3, 3, QColor(255,0,0), plot_data);
    ploter_win.drawPlot(plot);
    ploter_win.drawPlot(new Sin_plot(-2, 2, -3, 3, QColor(0,0,255)));
    */
	return a.exec();
}
