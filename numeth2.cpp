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
                 ((x_max() - x_min()) / (m_data.size() - 1));

        return m_data[id];
    }
private:
    std::vector<double> m_data;
};

class PieceLinPlot : public Plot
{
public:
    PieceLinPlot( double xMin, double xMax,
                  double yMin, double yMax,
                  QColor color,
                  std::map< double, double > const & grid )
        : Plot( xMin, xMax, yMin, yMax, color )
        , m_grid( grid )
    {}
    PieceLinPlot( double xMin, double xMax,
                  double yMin, double yMax,
                  QColor color,
                  std::vector< double > const & values )
        : Plot( xMin, xMax, yMin, yMax, color )
    {
        double step = ( ( xMax - xMin ) / ( values.size( ) - 1) );
        double x = xMin;
        for( auto val : values )
        {
            m_grid[ x ] = val;
            x += step;
        }
    }
    double get_val( double x )
    {
        auto nextNode = m_grid.upper_bound( x );
        double xLeft = std::prev(nextNode)->first;
        double yLeft = std::prev(nextNode)->second;
        double xRight = (nextNode)->first;
        double yRight = (nextNode)->second;
        if( nextNode == m_grid.end() )
        {
            return yLeft;
        }

        return ( yLeft +
                 ( x - xLeft ) / ( xRight - xLeft ) *
                 ( yRight - yLeft ) );
    }
private:
    std::map< double, double > m_grid;
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
    return cos(2*PI*x);
}

double sol_x_1(double x)
{
    return exp(x)*cos(2*PI*(x + 1));
}

double right(double x, double y)
{
    return exp(x*y)*(x*x + y*y - 8*PI*PI)*cos(2*PI*(x + y))
           - 4*PI*exp(x*y)*(x+y)*sin(2*PI*(x + y));
}

void expl(AnimPloter & ploter_win, std::ostream & fout)
{
    int M = 2056;
    int N_x = 64;
    int N_y = 64;
    double tau = 1./M;
    double h_x = 1./N_x;
    //fout << " " << h_x << std::endl;
    double h_y = 1./N_y;

    double eps = 0.001;
    double max = 0;
    double error = 0;
    int time = 0;

    std::vector<std::vector<double>> U_prev(N_x + 1, std::vector<double>(N_y + 1, 2));
    std::vector<std::vector<double>> U_new(N_x + 1, std::vector<double>(N_y + 1, 0));

    for(int j = 0; j <= N_y; ++j)
    {
        U_prev[0][j] = sol_0_y(j*h_y);
        U_prev[N_x][j] = sol_1_y(j*h_y);
        //U_prev[0][j] = 1;
        //U_prev[N_x][j] = 1;
    }
    for(int i = 0; i <= N_x; ++i)
    {
        U_prev[i][0] = sol_x_0(i*h_x);
        U_prev[i][N_y] = sol_x_1(i*h_x);
        //U_prev[i][0] = 1;
        //U_prev[i][N_y] = 1;
    }

    while(true)
    {
        max = 0;
        ++time;

        for(int j = 0; j <= N_y; ++j)
        {
            U_new[0][j] = sol_0_y(j*h_y);
            U_new[N_x][j] = sol_1_y(j*h_y);
            //U_new[0][j] = 1;
            //U_new[N_x][j] = 1;
        }
        for(int i = 0; i <= N_x; ++i)
        {
            U_new[i][0] = sol_x_0(i*h_x);
            U_new[i][N_y] = sol_x_1(i*h_x);
            //U_new[i][0] = 1;
            //U_new[i][N_y] = 1;
        }

        for(int i = 1; i < N_x; ++i)
        {
            for(int j = 1; j < N_y; ++j)
            {
                U_new[i][j] = tau/h_x/h_x*(U_prev[i-1][j] + U_prev[i+1][j])
                              + tau/h_y/h_y*(U_prev[i][j-1] + U_prev[i][j+1])
                              - tau*right(i*h_x, j*h_y)
                              + (1 - 2*tau/h_x/h_x -2*tau/h_y/h_y)*U_prev[i][j];
            }
        }

        for(int i = 1; i < N_x; ++i)
        {
            for(int j = 1; j < N_y; j++)
            {
                if(max < fabs((U_new[i][j] - U_prev[i][j])/tau))
                {
                    max = fabs((U_new[i][j] - U_prev[i][j])/tau);
                }
            }
        }
        std::cout << "max : " << max << std::endl;

        if(max < eps) break;

        U_prev = U_new;

    }

    for (int i = 0; i <= N_x; ++i)
    {
        for (int j = 0; j <= N_y; ++j)
        {
            if(fabs(U_new[i][j] - solution2(i*h_x, j*h_y)) > error)
            {
                error = fabs(U_new[i][j] - solution2(i*h_x, j*h_y));

            }
        }
    }
    fout << " " << error << std::endl;

    std::vector< Plot * > plot_anim;
    std::vector< Plot * > sol_plot_anim;
    for (int i = 0; i <= N_x; ++i)
    {
       /* for(auto el : U_new[i])
        {
            std::cout << el << '\t';
        }
        std::cout << std::endl;*/

        plot_anim.push_back(new PieceLinPlot(0, 1, -2, 2, QColor(255, 64, 64), U_new[i]));
        sol_plot_anim.push_back(new U_plot(0, 1, -2, 2, QColor(64, 64, 255), i*h_x));
    }
    ploter_win.drawAnimPlot(plot_anim);
    ploter_win.drawAnimPlot(sol_plot_anim);
    fout << "\n\n\n " << time << std::endl;
}
/*
void expl_appr(AnimPloter & ploter_win, std::ostream & fout_appr)
{
    int M = 100;
    int N_x_MAX = 1000;
    int N_y = 10;
    double tau = 1./M;
    //fout << " " << h_x << std::endl;
    double h_y = 1./N_y;

    double eps = 0.001;


    for (int N_x = 100; N_x <= N_x_MAX; N_x *= 2)
    {
        double h_x = 1./N_x;
        double max = 0;
        double error = 0;

        std::vector<std::vector<double>> U_prev(N_x + 1, std::vector<double>(N_y + 1, 0));
        std::vector<std::vector<double>> U_new(N_x + 1, std::vector<double>(N_y + 1, 0));

        for(int j = 0; j <= N_y; ++j)
        {
            U_prev[0][j] = sol_0_y(j*h_y);
            U_prev[N_x][j] = sol_1_y(j*h_y);
        }
        for(int i = 0; i <= N_x; ++i)
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
            for(int i = 0; i <= N_x; ++i)
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
                                  - tau*right(i*h_x, j*h_y) + (1 - 2*tau/h_x/h_x -2*tau/h_y/h_y)*U_prev[i][j];
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


        for (int i = 2; i <= N_x - 20; ++i)
        {
            for (int j = 1; j <= N_y - 1; ++j)
            {
                double new_err = fabs(U_new[i][j] - solution2(i*h_x, j*h_y));
                std::cout << U_new[i][j] << "\t"
                          << solution2(i*h_x, j*h_y) << "\t"
                          << new_err << "\t"
                          << i*h_x << "\t"
                          << j*h_y << std::endl;
                if(new_err > error)
                {
                    error = new_err;
                }

                }
            }
        }
        fout_appr << " " << error << std::endl;
    }

*/

void triangle(AnimPloter & ploter_win, std::ostream & fout)
{
    int M = 10240;
    int N_x = 32;
    int N_y = 32;
    double tau = 1./M;
    double h_x = 1./N_x;
    double h_y = 1./N_y;

    double eps = 0.001;
    double max = 0;
    double error = 0;
    int time = 0;

    std::vector<std::vector<double>> U_prev(N_x + 1, std::vector<double>(N_y + 1, 0));
    std::vector<std::vector<double>> U_new(N_x + 1, std::vector<double>(N_y + 1, 0));
    std::vector<std::vector<double>> Half_Xi(N_x + 1, std::vector<double>(N_y + 1, 0));
    std::vector<std::vector<double>> Xi(N_x + 1, std::vector<double>(N_y + 1, 0));

    for(int j = 0; j <= N_y; ++j)
    {
        U_prev[0][j] = sol_0_y(j*h_y);
        U_prev[N_x][j] = sol_1_y(j*h_y);
        //U_prev[0][j] = 1;
        //U_prev[N_x][j] = 1;
    }
    for(int i = 0; i <= N_x; ++i)
    {
        U_prev[i][0] = sol_x_0(i*h_x);
        U_prev[i][N_y] = sol_x_1(i*h_x);
        //U_prev[i][0] = 1;
        //U_prev[i][N_y] = 1;
    }

    while(true)
    {
        max = 0;
        ++time;
        for(int i = 0; i <= N_x; ++i)
        {
            for(int j = 0; j <= N_y; ++j)
            {
                Xi[i][j] = 0;
                Half_Xi[i][j] = 0;
            }
        }

        for(int j = 0; j <= N_y; ++j)
        {
            U_new[0][j] = sol_0_y(j*h_y);
            U_new[N_x][j] = sol_1_y(j*h_y);
            //U_new[0][j] = 1;
            //U_new[N_x][j] = 1;
        }
        for(int i = 0; i <= N_x; ++i)
        {
            U_new[i][0] = sol_x_0(i*h_x);
            U_new[i][N_y] = sol_x_1(i*h_x);
            //U_new[i][0] = 1;
            //U_new[i][N_y] = 1;
        }

        for (int i = 1; i < N_x; ++i)
        {
            for (int j = 1; j < N_y; ++j)
            {
                Half_Xi[i][j] = (tau/h_x/h_x*U_prev[i-1][j] + tau/h_x/h_x*U_prev[i+1][j]
                                + tau/h_y/h_y*U_prev[i][j-1] + tau/h_y/h_y*U_prev[i][j+1]
                                - 2*tau*(1/h_x/h_x + 1/h_y/h_y)*U_prev[i][j]
                                + tau/h_x/h_x * Half_Xi[i-1][j] + tau/h_y/h_y*Half_Xi[i][j-1]
                                - tau*right(i*h_x, j*h_y))
                                / (1 + tau/h_x/h_x + tau/h_y/h_y);
            }
        }

        for(int i = N_x - 1; i > 0; --i)
        {
            for(int j = N_y-1; j > 0; --j)
            {
                //Xi[i+1][j] = (Xi[i][j]*(1 + tau/h_x/h_x + tau/h_y/h_y)
                //             - tau/h_y/h_y*Xi[i][j+1] - Half_Xi[i][j])/(tau/h_x/h_x);
                //Xi[i][j+1] = (Xi[i][j]*(1 + tau/h_x/h_x + tau/h_y/h_y)
                //             - tau/h_x/h_x*Xi[i+1][j] - Half_Xi[i][j])/(tau/h_y/h_y);
                Xi[i][j] = (Half_Xi[i][j] + tau/h_x/h_x*Xi[i+1][j] + tau/h_y/h_y*Xi[i][j+1])/
                           (1 + tau/h_x/h_x + tau/h_y/h_y);
            }
        }

        /*
        for (int i = 1; i < N_x; ++i)
        {
            for (int j = 1; j < N_y; ++j)
            {
                Half_Xi[i+1][j] = (Half_Xi[i][j]*(1 + tau/h_x/h_x + tau/h_y/h_y) - tau/h_y/h_y*Half_Xi[i][j+1]
                                  - (tau/h_x/h_x*U_prev[i-1][j] + tau/h_x/h_x*U_prev[i+1][j]
                                  + tau/h_y/h_y*U_prev[i][j-1] + tau/h_y/h_y*U_prev[i][j+1]
                                  - 2*tau*(1/h_x/h_x + 1/h_y/h_y)*U_prev[i][j]))/(tau/h_y/h_y);
            }
        }

        for(int i = 1; i < N_x - 1; ++i)
        {
            for(int j = 1; j < N_y; ++j)
            {
                Xi[i][j] = Half_Xi[i][j] - tau/h_x/h_x*Xi[i-1][j] - tau/h_y/h_y*Xi[i][j-1];
            }
        }
*/
        for(int i = 1; i < N_x; ++i)
        {
            for(int j = 1; j < N_y; ++j)
            {
                U_new[i][j] = U_prev[i][j] + Xi[i][j];
            }
        }

        for(int i = 1; i < N_x; ++i)
        {
            for(int j = 1; j < N_y; j++)
            {
                if(max < fabs((U_new[i][j] - U_prev[i][j])/tau))
                {
                    max = fabs((U_new[i][j] - U_prev[i][j])/tau);
                }
            }
        }
        std::cout << "max : " << max << std::endl;

        if(max < eps) break;

        U_prev = U_new;
    }

    std::vector< Plot * > plot_anim;
    std::vector< Plot * > sol_plot_anim;
    for (int i = 0; i <= N_x; ++i)
    {
       /* for(auto el : U_new[i])
        {
            std::cout << el << '\t';
        }
        std::cout << std::endl;*/

        plot_anim.push_back(new PieceLinPlot(0, 1, -2, 2, QColor(255, 64, 64), U_new[i]));
        sol_plot_anim.push_back(new U_plot(0, 1, -2, 2, QColor(64, 64, 255), i*h_x));
    }
    ploter_win.drawAnimPlot(plot_anim);
    ploter_win.drawAnimPlot(sol_plot_anim);
    fout << "\n\n\n " << time << std::endl;

}
