#include "numeth.h"

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

std::vector< double > shuttle(std::vector< double > sub,
                              std::vector< double > main,
                              std::vector< double > super,
                              std::vector< double > right)
{
    int size = sub.size();
    std::vector< double > result(size, 0.);

    for (int i = 1; i < size; i++)
    {
        main[i] -= sub[i-1] * super[i-1] / main[i-1];
        right[i] -= sub[i-1] * right[i-1] / main[i-1];
    }

    result[size-1] = right[size-1] / main[size-1];
    for (int i = size - 2; i >= 0; i--)
        result[i] = (right[i] - super[i] * result[i+1]) / main[i];

    return result;
}

void task1(AnimPloter & ploter_win, std::ostream & fout)
{
    int N = 10;
    int M = 40;

    double alph = 1;

    double h = 1./(double)N;
    double tau = 1./(double)M;

    std::vector<double> first_row(N, 0.);

    std::vector< double > prev_row = first_row;
    std::vector< double > new_row(N, 0.);

    std::vector< Plot * > plot_anim;

    for(int t = 1; t < M; ++t)
    {
        // делаем график
        plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), prev_row));
        // печатаем в файл
        for(auto val : prev_row)
        {
            fout << val << " ";
        }
        fout << std::endl;

        new_row[0] = 1;

        std::vector< double > sub_coeff(N - 1, (alph / (2 * h)));
        std::vector< double > main_coeff(N, (2 * alph) / h);
        std::vector< double > super_coeff(N - 1, (1. / tau + (3 * alph) / (2 * h)));
        sub_coeff[N-2] = -1;
        main_coeff[N-1] = 1;

        std::vector< double > right_coeff(N, 0.);
        right_coeff[0] = (1. / tau - (1 - alph) / (2 * h)) * prev_row[2] +
                         2*(1 - alph) / h * prev_row[1] -
                         (1 - alph) / (2 * h) * prev_row[0] -
                         alph / (2 * h) * new_row[0];
        for (int j = 3; j <= N; ++j)
        {
            right_coeff[j-2] = (1. / tau - (1 - alph) / (2 * h)) * prev_row[j] +
                               2*(1 - alph) / h * prev_row[j-1] -
                               (1 - alph) / (2 * h) * prev_row[j-2];
        }
        right_coeff[N-1] = prev_row[N];

        std::vector< double > result = shuttle(sub_coeff, main_coeff, super_coeff, right_coeff);
        prev_row[0] = new_row[0];
        for(int j = 1; j <= N; ++j)
        {
            prev_row[j] = result[j-1];
        }
    }

    ploter_win.drawAnimPlot(plot_anim);
}
