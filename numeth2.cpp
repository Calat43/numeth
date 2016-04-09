#include "numeth2.h"

const double PI = 3.14159265;

double solution2(double x, double y)
{
    return exp(x*y)*cos(2.*PI*(x + y));
}

class U_plot : public Plot
{
public:
    U_plot(double y_min, double y_max,
             double u_min, double u_max,
             QColor color,
             double x)
        : Plot(y_min, y_max, u_min, u_max, color)
        , m_param(x)
    {}
    double get_val(double y)
    {
        return solution2(m_param, y);
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

double sol_0_y(double y)
{
    return cos(2*PI*y);
}

double sol_1_y(double y)
{
    return exp(y)*cos(2*PI*(y + 1));
}

double sol_x_0(double x)
{
    return cos(2*PI*(x));
}

double sol_x_1(double x)
{
    return exp(x)*cos(2*PI*(x + 1));
}

double right(double x, double y)
{
    return exp(x*y)*(x*x + y*y - 8*PI*PI)*cos(2*PI*(x + y))
           - 4*PI*exp(x*y)*sin(2*PI*(x + y));
}

void expl(AnimPloter & ploter_win)
{
    int M = 100;
    int N_x = 1000;
    int N_y = 1000;
    double tau = 1./M;
    double h_x = 1./N_x;
    double h_y = 1./N_y;

    double eps = 0.001;
    double max = 0;

    std::vector<std::vector<double>> U_prev(N_x + 1, std::vector<double>(N_y + 1, 0));
    std::vector<std::vector<double>> U_new(N_x + 1, std::vector<double>(N_y + 1, 0));

    for(int j = 0; j <= N_y; ++j)
    {
        U_prev[0][j] = sol_0_y(j*h_y);
        U_prev[N_x][j] = sol_1_y(j*h_y);
    }
    for(int i = 0; i <= N_y; ++i)
    {
        U_prev[i][0] = sol_x_0(i*h_x);
        U_prev[i][N_y] = sol_x_1(i*h_x);
    }

    while(true)
    {
        for(int j = 0; j <= N_y; ++j)
        {
            U_new[0][j] = sol_0_y(j*h_y);
            U_new[N_x][j] = sol_1_y(j*h_y);
        }
        for(int i = 0; i <= N_y; ++i)
        {
            U_new[i][0] = sol_x_0(i*h_x);
            U_new[i][N_y] = sol_x_1(i*h_x);
        }

        for(int i = 1; i < N_x; ++i)
        {
            for(int j = 1; j < N_y; ++j)
            {
                U_new[i][j] = tau/h_x/h_x*(U_prev[i-1][j] + U_prev[i+1][j])
                              + tau/h_y/h_y*(U_prev[i][j-1] + U_prev[i][j+1])
                              - tau*right(i*h_x, j*h_y) + (M - 2*tau/h_x/h_x -2*tau/h_y/h_y)*U_prev[i][j];
            }
        }

        for(int i = 0; i <= N_x; ++i)
        {
            for(int j = 0; j <= N_y; j++)
            {
               max = fabs((U_new[i][j] - U_prev[i][j])/tau);
            }
        }

        if(max < eps) break;

        U_prev = U_new;

    }

    std::vector< Plot * > plot_anim;
    std::vector< Plot * > sol_plot_anim;
    for (int i = 0; i <= N_x; ++i)
    {
        plot_anim.push_back(new Array_plot(0, 1, -2, 2, QColor(255, 64, 64), U_new[i]));
        sol_plot_anim.push_back(new U_plot(0, 1, -2, 2, QColor(64, 64, 255), i*h_x));
    }
    ploter_win.drawAnimPlot(plot_anim);
    ploter_win.drawAnimPlot(sol_plot_anim);
}
