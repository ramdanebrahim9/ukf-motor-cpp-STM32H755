#pragma once

static inline double UKF_interp1(const double* x_table,
                                 const double* y_table,
                                 int len,
                                 double x_query)
{
    if (x_query <= x_table[0])
        return y_table[0];
    if (x_query >= x_table[len - 1])
        return y_table[len - 1];

    for (int i = 0; i < len - 1; i++)
    {
        if (x_query <= x_table[i + 1])
        {
            double t = (x_query - x_table[i]) / (x_table[i + 1] - x_table[i]);
            return y_table[i] + t * (y_table[i + 1] - y_table[i]);
        }
    }
    return y_table[len - 1];
}