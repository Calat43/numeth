#include "numeth.h"

double solution(double x, double t)
{
    if (t < x) return 0;
    return 1;
}

class Sin_plot : public Plot
{
public:
    Sin_plot(double x_min, double x_max,
             double y_min, double y_max,
             QColor color,
             double t)
        : Plot(x_min, x_max, y_min, y_max, color)
        , m_param(t)
    {}
    double get_val(double x)
    {
        return sin(m_param - x);
    }
private:
    double m_param;
};

class U_plot : public Plot
{
public:
    U_plot(double x_min, double x_max,
             double y_min, double y_max,
             QColor color,
             double t)
        : Plot(x_min, x_max, y_min, y_max, color)
        , m_param(t)
    {}
    double get_val(double x)
    {
        return solution(x, m_param);
    }
private:
    double m_param;
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
    int const N_MAX = 1000;

    int N = 100;
    int M = 4000;

    double alph = 0;


    std::vector<double> first_row(N, 0.);

    std::vector< double > prev_row = first_row;
    std::vector< double > new_row(N, 0.);

    std::vector< Plot * > plot_anim;

    for (int N = 2; N <= N_MAX; N*=2) {
        double h = 1./(double)N;
        double tau = 1./(double)M;

        double err = 0;
        //бегущий счет для альфа = 0
        for(int t = 0; t <= M; ++t)
        {
            prev_row.resize(N);
            new_row.resize(N);
            //plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), prev_row));
            new_row[0] = 1;
            new_row[1] = 1;
            /*for(auto val : prev_row)
            {
                fout << val << " ";
            }
            fout << std::endl;*/
            for(int j = 2; j <= N; ++j)
            {
                new_row[j] = (1. - 1.5*tau*(double)N)*prev_row[j] +
                             (2.*tau*(double)N)*prev_row[j-1] - (0.5*tau*(double)N)*prev_row[j-2];
             }

            for(int j = 0; j <= N; ++j)
            {
                prev_row[j] = new_row[j];
            }
          //  if (t == 3*M/5) {
            for (int j = 0; j <= N; ++j)
            {
                double new_err = fabs(new_row[j] - solution(j*h, t*tau));
                if (err < new_err) err = new_err;
            }
          //  }
        }
        fout << N << " " << err << std::endl;
    }
/*
    //прогонка для общего случая
    for(int t = 1; t < M; ++t)
    {
        // делаем график
        //plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), prev_row));
        plot_anim.push_back(new U_plot(-1, 2, -1, 1, QColor(255,0,0), t*tau));
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
*/

    ploter_win.drawAnimPlot(plot_anim);
}

void sin_task1(AnimPloter & ploter_win, std::ostream & fout)
{
    int const N_MAX = 1000;

    int N = 100;
    int M = 4000;

    double alph = 0;

    double h = 1./(double)N;
    double tau = 1./(double)M;

    std::vector<double> first_row(N, 0.);

    for (int j = 0; j < N; ++j)
    {
        first_row[j] = -sin(j*h);
    }

    std::vector< double > prev_row = first_row;
    std::vector< double > new_row(N, 0.);

    std::vector< Plot * > plot_anim;

    //бегущий счет для альфа = 0
    for(int t = 0; t <= M; ++t)
    {
        prev_row.resize(N);
        new_row.resize(N);
        //plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), prev_row));
        plot_anim.push_back(new Sin_plot(0, 1, -2, 2, QColor(255, 64, 64), t*tau));
        new_row[0] = sin(t*tau);
        new_row[1] = sin(t*tau - h);

        for(int j = 2; j <= N; ++j)
        {
            new_row[j] = (1. - 1.5*tau*(double)N)*prev_row[j] +
                         (2.*tau*(double)N)*prev_row[j-1] - (0.5*tau*(double)N)*prev_row[j-2];
        }

        for(int j = 0; j <= N; ++j)
        {
            prev_row[j] = new_row[j];
        }
    }

    /* for (int N = 2; N <= N_MAX; N*=2) {
        double h = 1./(double)N;
        double tau = 1./(double)M;

        double err = 0;
        //бегущий счет для альфа = 0
        for(int t = 0; t <= M; ++t)
        {
            prev_row.resize(N);
            new_row.resize(N);
            //plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), prev_row));
            new_row[0] = 1;
            new_row[1] = 1;

            for(int j = 2; j <= N; ++j)
            {
                new_row[j] = (1. - 1.5*tau*(double)N)*prev_row[j] +
                             (2.*tau*(double)N)*prev_row[j-1] - (0.5*tau*(double)N)*prev_row[j-2];
             }

            for(int j = 0; j <= N; ++j)
            {
                prev_row[j] = new_row[j];
            }
          //  if (t == 3*M/5) {
            for (int j = 0; j <= N; ++j)
            {
                double new_err = fabs(new_row[j] - solution(j*h, t*tau));
                if (err < new_err) err = new_err;
            }
          //  }
        }
        fout << N << " " << err << std::endl;
    }
   */

    ploter_win.drawAnimPlot(plot_anim);
}
