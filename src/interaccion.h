#include <stdio.h>
#include <stdlib.h>


float min_diff(float x, float y, float L);
float r_squared(float* x_1, float* x_2, float L);
float potential(float r_squared, float sigma);
float get_force_from_table(float r2, float r_c, float *table_f, float *table_r, int length);
int update_forces(float *f, float *x, int N, float L, float r_c, float *table_f, float *table_r2, int length);
