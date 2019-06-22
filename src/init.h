#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592


float sample_boltzmann(float temperature);
int amount_of_particles(float rho, int size);
float initial_positions(int size, float* x, int n);
int initial_velocities(float* v, int n, float temperature);
int fill_forces_table(float *table_r, float *table_r2, float *table_f, float *table_v, float r_c, int length);
