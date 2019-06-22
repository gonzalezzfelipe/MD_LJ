#include <stdio.h>
#include <stdlib.h>


float min_diff(float x, float y, int periodicity);
float r_squared(float* x_1, float* x_2, int size);
float potential(float r_squared, float sigma);
float get_force_from_table(float r2, float r_c, float *table_f, float *table_r, int length);
int update_forces(float *f, float *x, int N, int L, float r_c, float *table_f, float *table_r2, int length);
