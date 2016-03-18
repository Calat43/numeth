#include "numeth.h"

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

void task1(std::vector< double > first_row)
{
    int N = first_row.size(); // better get as arg
    int M = 8; // same

    std::vector< double > prev_row = first_row;
    std::vector< double > new_row(N, 0.);

    for(int t = 1; t < M; ++t)
    {
        new_row[0] = 1;

        double h = 1;
        double tau = 1;
        double alph = 1;

        std::vector< double > sub_coeff(N - 1, (alph / (2 * h)));
        std::vector< double > main_coeff(N, (2 * alph) / h);
        std::vector< double > super_coeff(N - 1, (1. / tau + (3 * alph) / (2 * h)));
        sub_coeff[N-2] = 1;
        main_coeff[N-1] = 1;

        std::vector< double > right_coeff(N, 0.);
        rigth_coeff[0] = (1. / tau - (1 - alph) / (2 * h)) * prev_row[2] +
                         2*(1 - alph) / h * prev_row[1] -
                         (1 - alph) / (2 * h) * prev_row[0] -
                         alph / (2 * h) * new_row[0];
        for (int j = 3; j <= N; ++j)
        {
            rigth_coeff[j-2] = (1. / tau - (1 - alph) / (2 * h)) * prev_row[j] +
                               2*(1 - alph) / h * prev_row[j-1] -
                               (1 - alph) / (2 * h) * prev_row[j-2];
        }
        rigth_coeff[N-1] = 0;

        std::vector< double > result = shuttle(sub_coeff, main_coeff, super_coeff, right_coeff);
        prew_row[0] = new_row[0];
        for(int j = 1; j <= N; ++j)
        {
            prew_row[j] = result[j-1];
        }
    }






}
