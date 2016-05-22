#include "numeth3.h"
/*
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
*/
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

double w0_1(double x, double x0, double ro_l, double ro_r)
{
    if(x <= x0) {return ro_l;}
    return ro_r;
}
double w0_2(double x, double x0, double ro_l, double ro_r, double u_l, double u_r)
{
    if(x <= x0)
    {

        return ro_l*u_l;
    }
    return ro_r*u_r;
}
double w0_3(double x, double x0, double ro_l, double ro_r, double u_l, double u_r, double p_l, double p_r, double y)
{
    if(x <= x0)
    {
        return p_l/(y - 1) + ro_l/2*u_l*u_l;
    }
    return p_r/(y - 1) + ro_r/2*u_r*u_r;
}

double e( double p, double ro, double u, double y )
{
    return p / (y - 1) + ro / 2 * u*u;
}
double p_from_e( double e, double ro, double u, double y )
{
    return (e - ro / 2 * u * u ) * (y - 1);
}

double f_1(  double w_1,
             double w_2,
             double w_3 )
{
    return w_2;
}
double f_2(  double w_1,
             double w_2,
             double w_3, double y )
{
    return w_2 * w_2 * p_from_e( w_3, w_1, w_2 / w_1, y );
}
double f_3(  double w_1,
             double w_2,
             double w_3, double y )
{
    double p = p_from_e( w_3, w_1, w_2 / w_1, y );
    return ( w_3 + p ) * w_2 / w_1;
}

double f_astr_1( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y )
{
    double f = f_1( w_1[j], w_2[j], w_3[j] );
    return f - D * h * ( w_1[j] - w_1[j - 1] );
}
double f_astr_2( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y )
{
    double f = f_2( w_1[j], w_2[j], w_3[j], y );
    return f - D * h * ( w_2[j] - w_2[j - 1] );
}
double f_astr_3( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y )
{
    double f = f_3( w_1[j], w_2[j], w_3[j], y );
    return f - D * h * ( w_3[j] - w_3[j - 1] );
}

double w_tild_1( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    return w_1[j] - tau / h *
            ( f_astr_1( w_1, w_2, w_3, j+1, D, h, y ) -
              f_astr_1( w_1, w_2, w_3, j, D, h, y ) );
}
double w_tild_2( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    return w_2[j] - tau / h*
            ( f_astr_2( w_1, w_2, w_3, j+1, D, h, y ) -
              f_astr_2( w_1, w_2, w_3, j, D, h, y ) );
}
double w_tild_3( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    return w_3[j] - tau / h *
            ( f_astr_3( w_1, w_2, w_3, j+1, D, h, y ) -
              f_astr_3( w_1, w_2, w_3, j, D, h, y ) );
}

double w_next_1( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double f = f_1( w_t_1, w_t_2, w_t_3 );
    double fjm1 = f_1( w_t_1jm1, w_t_2jm1, w_t_3jm1 );
    return (w_1[j] + w_t_1) / 2 - tau/2/h*( f - fjm1 );
}
double w_next_2( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double f    = f_2( w_t_1, w_t_2, w_t_3, y );
    double fjm1 = f_2( w_t_1jm1, w_t_2jm1, w_t_3jm1, y );
    return (w_2[j] + w_t_2) / 2 - tau/2/h*( f - fjm1 );
}
double w_next_3( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double f = f_3( w_t_1, w_t_2, w_t_3, y );
    double fjm1 = f_3( w_t_1jm1, w_t_2jm1, w_t_3jm1, y );
    return (w_3[j] + w_t_3) / 2 - tau/2/h*( f - fjm1 );
}

void maccormack( AnimPloter & ploter_win,
                 std::ostream & fout )
{
    //показатель политропы
    double y = 1.4;
    int M = 80;
    int N = 10;
    double T = 1;
    double D = 5;
    //начальные данные, p_l > p_r
    double ro_l = 8;
    double p_l = 10/y;
    double u_l = 0;
    double ro_r = 1;
    double p_r = 1/y;
    double u_r = 0;
    double xl = 0;
    double xr = 1;
    double x0 = 0.5;

    double h = (xr - xl)/(M + 1);
    double tau = T/(N + 1);

    std::vector<std::vector<double>> ro(M+1, std::vector<double>(N+1, 0));
    std::vector<std::vector<double>> u(M+1, std::vector<double>(N+1, 0));
    std::vector<std::vector<double>> p(M+1, std::vector<double>(N+1, 0));

    std::vector< double > w_prev_1(N+1, 0);
    std::vector< double > w_prev_2(N+1, 0);
    std::vector< double > w_prev_3(N+1, 0);
    for( int j = 0; j <= N; ++j)
    {
        w_prev_1[j] = w0_1( j*h, x0, ro_l, ro_r );
        w_prev_2[j] = w0_2( j*h, x0, ro_l, ro_r, u_l, u_r );
        w_prev_3[j] = w0_3( j*h, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );
    }
    for( int t = 0; t <= M; ++t )
    {
        std::vector< double > w_new_1(N+1, 0);
        std::vector< double > w_new_2(N+1, 0);
        std::vector< double > w_new_3(N+1, 0);

        w_new_1[0] = w0_1( 0, x0, ro_l, ro_r );
        w_new_2[0] = w0_2( 0, x0, ro_l, ro_r, u_l, u_r );
        w_new_3[0] = w0_3( 0, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );
        for( int j = 1; j <= N - 1; ++j)
        {
            w_new_1[j] = w_next_1( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
            w_new_2[j] = w_next_2( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
            w_new_3[j] = w_next_3( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
        }
        w_new_1[M] = w0_1( (M)*h, x0, ro_l, ro_r );
        w_new_2[M] = w0_2( (M)*h, x0, ro_l, ro_r, u_l, u_r );
        w_new_3[M] = w0_3( (M)*h, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );

        for( int j = 0; j <= N; ++j)
        {
            ro[t][j] = w_new_1[j];
            u[t][j] = w_new_2[j] / w_new_1[j];
            p[t][j] = p_from_e( w_new_3[j], ro[t][j], u[t][j], y );
        }
        std::swap( w_prev_1, w_new_1 );
        std::swap( w_prev_2, w_new_2 );
        std::swap( w_prev_3, w_new_3 );
    }

    std::vector< Plot * > ro_plot;
    std::vector< Plot * > sol_ro_plot;
    std::vector< Plot * > u_plot;
    std::vector< Plot * > sol_u_plot;
    std::vector< Plot * > p_plot;
    std::vector< Plot * > sol_p_plot;
    for (int t = 0; t <= M; ++t)
    {
        for( int j = 0; j <= N; ++j)
        {
            fout << ro[t][j] << " ";
        }
        fout << std::endl;
        ro_plot.push_back(new PieceLinPlot(0, 1, -2, 2, QColor(255, 64, 64), ro[t]));
        u_plot.push_back(new PieceLinPlot(0, 1, -2, 2,  QColor(64, 255, 64), u[t]));
        p_plot.push_back(new PieceLinPlot(0, 1, -2, 2,  QColor(64, 64, 255), p[t]));
        //sol_plot_anim.push_back(new U_plot(0, 1, -2, 2, QColor(64, 64, 255), i*h_x));
    }
    ploter_win.drawAnimPlot(ro_plot);
    ploter_win.drawAnimPlot(u_plot);
    ploter_win.drawAnimPlot(p_plot);
}
