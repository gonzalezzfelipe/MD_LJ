#include <stdio.h>
#include <stdlib.h>


// Termalizacion
float verlet_coeff(float *x, float L, int N);
float h_boltzmann(float *v, int N);
float g_r(float *x, int N);

// Termodinamica
float potential_energy(float *x, float *v, float *table_r2, float *table_v, int N, float L, float r_c, int lenght);
float temperature(float *v, int N);
float pressure(float rho, float T, float L, float *table_f, float *table_r, float r_c, int lenght, float *x, int N);

// Escribir
int write_log(
    int timestep,
    char* filename,
    float rho,
    int N,
    float L,
    float r_c,
    float *table_r2,
    float *table_v,
    float *table_f,
    int lenght,
    float *x,
    float *v,
    float *f);
