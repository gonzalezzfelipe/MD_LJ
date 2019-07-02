#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592


double sample_boltzmann(double temperature);
int amount_of_particles(double rho, double L);
double initial_positions(double L, double* x, int n);
int initial_velocities(double* v, int n, double temperature);
int fill_forces_table(double *table_r, double *table_r2, double *table_f, double *table_v, double r_c, int length);

int rescaling(
    double T,
    double relative_error,
    int termalization,
    double* x,
    double* v,
    double* f,
    int N,
    double dt,
    double L,
    double r_c,
    double *table_f,
    double *table_r2,
    int length);
