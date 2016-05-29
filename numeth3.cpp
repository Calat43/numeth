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


double p3_func(double p3, double y, double p1, double ro1, double p2, double ro2, double u2)
{
    double left_part = pow(
                    (2/(y-1)*sqrt(y*pow(p1, 1/y)/ro1)*
                     (pow(p1, (y-1)/2/y) - pow(p3, (y-1)/2/y))
                    ) - u2,
                   2);
    double right_part = (p3 - p2) / ro2 *
                (1 - ((y-1)*p3 + (y+1)*p2) / ((y+1)*p3 + (y-1)*p2));

    return left_part - right_part;
}

double p3_bissection(double y, double p1, double ro1, double p2, double ro2, double u2, double eps)
{
    double left = p2;
    double right = p1;
    double dp = right - left;
    while(dp > eps)
    {
        dp = dp / 2;
        double middle = left + dp;
        if( p3_func(left, y, p1, ro1, p2, ro2, u2) *
            p3_func(middle, y, p1, ro1, p2, ro2, u2) < 0 )
        {
            right = middle;
        }
        else
        {
            left = middle;
        }
        dp = right - left;
    }
    return (left + right) / 2;
}

double get_u3( double p1, double ro1, double p2, double ro2, double y, double p3 )
{
    return 2/(y-1)*sqrt(y*pow(p1, 1/y)/ro1)*
            (pow(p1, (y-1)/2/y) - pow(p3, (y-1)/2/y));
}

double get_ro3( double p2, double p3, double ro2, double y )
{
    return ro2*(((y+1)*p3 + (y-1)*p2) / ((y-1)*p3 + (y+1)*p2));
}

double get_D( double ro2, double ro3, double u3 )
{
    return u3*ro3/(ro3 - ro2);
}

double exact_dencity( double x, double t, double ro1, double ro2, double p1, double p2, double y, double u2, double eps )
{
    double p3 = p3_bissection( y, p1, ro1, p2, ro2, u2, eps);
    double u3 = get_u3( p1, ro1, p2, ro2, y, p3);
    double ro3 = get_ro3( p2, p3, ro2, y );
    double D = get_D( ro2, ro3, u3 );
    double c1 = sqrt(y*p1/ro1);
    double c3_ = sqrt(y*p3/ro3);

    double x1 = -t * c1;
    double x2 = t*(u3 - c3_);
    double x3 = u3 * t;
    double x4 = D * t;

    if( x <= x1 )
    {
        return ro1;
    }
    if( x <= x2 )
    {
        double c = 2/(y + 1)*(c1 - (y-1)/2 * x/t );
        return pow( c*c/y * pow(ro1, y)/p1, 1/(y-1) );
    }
    if( x <= x3 )
    {
        return ro1 * pow( p3/p1, 1/y );
    }
    if( x <= x4 )
    {
        return ro3;
    }
    return ro2;
}

double exact_velocity( double x, double t, double ro1, double ro2, double p1, double p2, double y, double u1, double u2, double eps )
{
    double p3 = p3_bissection( y, p1, ro1, p2, ro2, u2, eps);
    double u3 = get_u3( p1, ro1, p2, ro2, y, p3);
    double ro3 = get_ro3( p2, p3, ro2, y );
    double D = get_D( ro2, ro3, u3 );
    double c1 = sqrt(y*p1/ro1);
    double c3_ = sqrt(y*p3/ro3);

    double x1 = -t * c1;
    double x2 = t*(u3 - c3_);
    double x3 = u3 * t;
    double x4 = D * t;

    if( x <= x1 )
    {
        return u1;
    }
    if( x <= x2 )
    {
        return 2 / (y + 1) * ( c1 + x/t );
    }
    if( x <= x3 )
    {
        return u3;
    }
    if( x <= x4 )
    {
        return u3;
    }
    return u2;
}

double exact_pressure( double x, double t, double ro1, double ro2, double p1, double p2, double y, double u2, double eps )
{
    double p3 = p3_bissection( y, p1, ro1, p2, ro2, u2, eps);
    double u3 = get_u3( p1, ro1, p2, ro2, y, p3);
    double ro3 = get_ro3( p2, p3, ro2, y );
    double D = get_D( ro2, ro3, u3 );
    double c1 = sqrt(y*p1/ro1);
    double c3_ = sqrt(y*p3/ro3);

    double x1 = -t * c1;
    double x2 = t*(u3 - c3_);
    double x3 = u3 * t;
    double x4 = D * t;

    if( x <= x1 )
    {
        return p1;
    }
    if( x <= x2 )
    {
        double c = 2/(y + 1)*(c1 - (y-1)/2 * x/t );
        double ro =  pow( c*c/y * pow(ro1, y)/p1, 1/(y-1) );
        return c*c/y*ro;
    }
    if( x <= x3 )
    {
        return p3;
    }
    if( x <= x4 )
    {
        return p3;
    }
    return p2;
}

class Dencity_Plot : public Plot
{
public:
    Dencity_Plot(double y_min, double y_max,
             double x_min, double x_max,
             QColor color,
             double t,
             double ro1, double ro2, double p1, double p2, double y, double u2, double eps )
        : Plot(y_min, y_max, x_min, x_max, color)
        , m_t(t)
        , m_ro1( ro1 )
        , m_ro2( ro2 )
        , m_p1( p1 )
        , m_p2( p2 )
        , m_y( y )
        , m_u2( u2 )
        , m_eps( eps )
    { }
    double get_val(double x)
    {
        return exact_dencity( x, m_t, m_ro1, m_ro2, m_p1, m_p2, m_y, m_u2, m_eps );
    }
private:
    double m_t;
    double m_ro1;
    double m_ro2;
    double m_p1;
    double m_p2;
    double m_y;
    double m_u2;
    double m_eps;
};

class Velocity_Plot : public Plot
{
public:
    Velocity_Plot(double y_min, double y_max,
             double x_min, double x_max,
             QColor color,
             double t,
             double ro1, double ro2, double p1, double p2, double y, double u1, double u2, double eps )
        : Plot(y_min, y_max, x_min, x_max, color)
        , m_t(t)
        , m_ro1( ro1 )
        , m_ro2( ro2 )
        , m_p1( p1 )
        , m_p2( p2 )
        , m_y( y )
        , m_u1( u1 )
        , m_u2( u2 )
        , m_eps( eps )
    { }
    double get_val(double x)
    {
        return exact_velocity( x, m_t, m_ro1, m_ro2, m_p1, m_p2, m_y, m_u1, m_u2, m_eps );
    }
private:
    double m_t;
    double m_ro1;
    double m_ro2;
    double m_p1;
    double m_p2;
    double m_y;
    double m_u1;
    double m_u2;
    double m_eps;
};

class Pressure_Plot : public Plot
{
public:
    Pressure_Plot(double y_min, double y_max,
             double x_min, double x_max,
             QColor color,
             double t,
             double ro1, double ro2, double p1, double p2, double y, double u2, double eps )
        : Plot(y_min, y_max, x_min, x_max, color)
        , m_t(t)
        , m_ro1( ro1 )
        , m_ro2( ro2 )
        , m_p1( p1 )
        , m_p2( p2 )
        , m_y( y )
        , m_u2( u2 )
        , m_eps( eps )
    { }
    double get_val(double x)
    {
        return exact_pressure( x, m_t, m_ro1, m_ro2, m_p1, m_p2, m_y, m_u2, m_eps );
    }
private:
    double m_t;
    double m_ro1;
    double m_ro2;
    double m_p1;
    double m_p2;
    double m_y;
    double m_u2;
    double m_eps;
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
    return w_2 * w_2 / w_1 + p_from_e( w_3, w_1, w_2 / w_1, y );
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

double f_tild_astr_1( double w_1, double w_2, double w_3,
                      double w_1p1, double w_2p1, double w_3p1,
                      double D, double h, double y )
{
    double f = f_1( w_1, w_2, w_3 );
    return f - D * h * ( w_1p1 - w_1 );
}
double f_tild_astr_2( double w_1, double w_2, double w_3,
                      double w_1p1, double w_2p1, double w_3p1,
                      double D, double h, double y )
{
    double f = f_2( w_1, w_2, w_3, y );
    return f - D * h * ( w_2p1 - w_2 );
}
double f_tild_astr_3( double w_1, double w_2, double w_3,
                      double w_1p1, double w_2p1, double w_3p1,
                      double D, double h, double y )
{
    double f = f_3( w_1, w_2, w_3, y );
    return f - D * h * ( w_3p1 - w_3 );
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
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_1jp1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2jp1 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3jp1 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double f = f_tild_astr_1( w_t_1, w_t_2, w_t_3, w_t_1jp1, w_t_2jp1, w_t_3jp1, D, h, y );
    double fjm1 = f_tild_astr_1( w_t_1jm1, w_t_2jm1, w_t_3jm1, w_t_1, w_t_2, w_t_3, D, h, y );
    return (w_1[j] + w_t_1) / 2 - tau/2/h*( f - fjm1 );
}
double w_next_2( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_1jp1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2jp1 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3jp1 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double f = f_tild_astr_2( w_t_1, w_t_2, w_t_3, w_t_1jp1, w_t_2jp1, w_t_3jp1, D, h, y );
    double fjm1 = f_tild_astr_2( w_t_1jm1, w_t_2jm1, w_t_3jm1, w_t_1, w_t_2, w_t_3, D, h, y );
    return (w_2[j] + w_t_2) / 2 - tau/2/h*( f - fjm1 );
}
double w_next_3( std::vector< double > w_1,
                 std::vector< double > w_2,
                 std::vector< double > w_3, int j,
                 double D, double h, double y, double tau )
{
    double w_t_1 = w_tild_1( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_2 = w_tild_2( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_3 = w_tild_3( w_1, w_2, w_3, j, D, h, y, tau );
    double w_t_1jm1 = w_tild_1( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_2jm1 = w_tild_2( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_3jm1 = w_tild_3( w_1, w_2, w_3, j-1, D, h, y, tau );
    double w_t_1jp1 = w_tild_1( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_2jp1 = w_tild_2( w_1, w_2, w_3, j+1, D, h, y, tau );
    double w_t_3jp1 = w_tild_3( w_1, w_2, w_3, j+1, D, h, y, tau );
    double f = f_tild_astr_3( w_t_1, w_t_2, w_t_3, w_t_1jp1, w_t_2jp1, w_t_3jp1, D, h, y );
    double fjm1 = f_tild_astr_3( w_t_1jm1, w_t_2jm1, w_t_3jm1, w_t_1, w_t_2, w_t_3, D, h, y );
    return (w_3[j] + w_t_3) / 2 - tau/2/h*( f - fjm1 );
}

void maccormack( AnimPloter & ro_ploter_win,
                 AnimPloter & u_ploter_win,
                 AnimPloter & p_ploter_win,
                 std::ostream & fout )
{
    double y = 1.4;
    int M = 500;
    int N = 200;
    double T = 2.5;
    double D = 4;
    //начальные данные, p_l > p_r
    double ro_l = 1;
    double p_l = 1;
    double u_l = 0;
    double ro_r = 0.125;
    double p_r = 0.1;
    double u_r = 0;
    double xl = -4;
    double xr = 6;
    double x0 = 0;

    double eps = 0.01;

    double h = (xr - xl)/(M + 1);
    double tau = T/(N + 1);

    std::vector<std::vector<double>> ro(N+1, std::vector<double>(M+1, 0));
    std::vector<std::vector<double>> u(N+1, std::vector<double>(M+1, 0));
    std::vector<std::vector<double>> p(N+1, std::vector<double>(M+1, 0));

    std::vector< double > w_prev_1(M+1, 0);
    std::vector< double > w_prev_2(M+1, 0);
    std::vector< double > w_prev_3(M+1, 0);
    for( int j = 0; j <= M; ++j)
    {
        w_prev_1[j] = w0_1( xl + j*h, x0, ro_l, ro_r );
        w_prev_2[j] = w0_2( xl + j*h, x0, ro_l, ro_r, u_l, u_r );
        w_prev_3[j] = w0_3( xl + j*h, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );
    }
    for( int t = 0; t <= N; ++t )
    {
        for( int j = 0; j <= M; ++j)
        {
            ro[t][j] = w_prev_1[j];
            u[t][j] = w_prev_2[j] / w_prev_1[j];
            p[t][j] = p_from_e( w_prev_3[j], ro[t][j], u[t][j], y );
        }

        std::vector< double > w_new_1(M+1, 0);
        std::vector< double > w_new_2(M+1, 0);
        std::vector< double > w_new_3(M+1, 0);

        w_new_1[0] = w0_1( xl, x0, ro_l, ro_r );
        w_new_2[0] = w0_2( xl, x0, ro_l, ro_r, u_l, u_r );
        w_new_3[0] = w0_3( xl, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );
        for( int j = 1; j <= M - 1; ++j)
        {
            w_new_1[j] = w_next_1( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
            w_new_2[j] = w_next_2( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
            w_new_3[j] = w_next_3( w_prev_1, w_prev_2, w_prev_3, j, D, h, y, tau );
        }
        w_new_1[M] = w0_1( xr, x0, ro_l, ro_r );
        w_new_2[M] = w0_2( xr, x0, ro_l, ro_r, u_l, u_r );
        w_new_3[M] = w0_3( xr, x0, ro_l, ro_r, u_l, u_r, p_l, p_r, y );
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
    for (int t = 0; t <= N; ++t)
    {
        for( int j = 0; j <= M; ++j)
        {
            fout << u[t][j] << " ";
        }
        fout << std::endl;
        ro_plot.push_back(new PieceLinPlot(-4, 6, 0, 1, QColor(255, 64, 64), ro[t]));
        u_plot.push_back(new PieceLinPlot(-4, 6, -1, 10,  QColor(255, 64, 64), u[t]));
        p_plot.push_back(new PieceLinPlot(-4, 6, 0, 1, QColor(255, 64, 64), p[t]));

        sol_ro_plot.push_back(new Dencity_Plot(-4, 6, 0, 1, QColor(64, 64, 255), t*tau, ro_l, ro_r, p_l, p_r, y, u_r, eps));
        sol_u_plot.push_back(new Velocity_Plot(-4, 6, 0, 1, QColor(64, 64, 255), t*tau, ro_l, ro_r, p_l, p_r, y, u_l, u_r, eps));
        sol_p_plot.push_back(new Pressure_Plot(-4, 6, 0, 1, QColor(64, 64, 255), t*tau, ro_l, ro_r, p_l, p_r, y, u_r, eps));
    }
    ro_ploter_win.drawAnimPlot(ro_plot);
    u_ploter_win.drawAnimPlot(u_plot);
    p_ploter_win.drawAnimPlot(p_plot);

    ro_ploter_win.drawAnimPlot(sol_ro_plot);
    u_ploter_win.drawAnimPlot(sol_u_plot);
    p_ploter_win.drawAnimPlot(sol_p_plot);
}
