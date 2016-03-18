#include "ploter.h"
#include <QApplication>
#include <QtMath>

class Sin_plot : public Plot
{
public:
	Sin_plot(double x_min, double x_max,
			 double y_min, double y_max,
			 QColor color)
		: Plot(x_min, x_max, y_min, y_max, color) {}
	double get_val(double x)
	{
		return qSin(1./x) / x;
	}
};

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Ploter ploter_win(-5, 5, -3, 3, QColor(0, 0, 0));
	ploter_win.show();
	Sin_plot * sin_plot =
			new Sin_plot(-2, 2, -2, 2,
						 QColor(255, 64, 64));
	ploter_win.drawPlot(sin_plot);

	return a.exec();
}
