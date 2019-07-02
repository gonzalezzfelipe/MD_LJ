#include <stdio.h>
#include <stdlib.h>


double min_diff(double x, double y, double L);
double r_squared(double* x_1, double* x_2, double L);
double potential(double r_squared, double sigma);
double get_force_from_table(double r2, double r_c, double *table_f, double *table_r, int length);
int update_forces(double *f, double *x, int N, double L, double r_c, double *table_f, double *table_r2, int length);
