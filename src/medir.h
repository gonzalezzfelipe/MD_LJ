#include <stdio.h>
#include <stdlib.h>


// Termalizacion
float verlet_coeff(float *x, int L, int N);
float h_boltzmann(float *v, int N);
float g_r(float *x, int N);

// Termodinamica
float potential_energy(float *x, float *v, float *table_r2, float *table_v, int N, int L, float r_c);
float temperature(float *v, int N);
float pressure(float rho, float T, int L, float *f, float *r);
