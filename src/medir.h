#include <stdio.h>
#include <stdlib.h>


// Termalizacion
double verlet_coeff(double *x, double L, int N);
double h_boltzmann(double *v, int N);
double g_r(double *x, int N);

// Termodinamica
double potential_energy(double *x, double *v, double *table_r2, double *table_v, int N, double L, double r_c, int lenght);
double temperature(double *v, int N);
double kinetic(double *v, int N);
double pressure(double rho, double T, double L, double *table_f, double *table_r, double r_c, int lenght, double *x, int N);

// Escribir
int write_log(
  int timestep,
  double time,
  char* filename,
  double rho,
  int N,
  double L,
  double r_c,
  double *table_r2,
  double *table_v,
  double *table_f,
  int lenght,
  double *x,
  double *v,
  double *f);
