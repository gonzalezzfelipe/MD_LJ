#include <stdio.h>
#include <stdlib.h>
#include "objects.h"


// Termalizacion
double verlet_coeff(Particles parts, double L);
double h_boltzmann(Particles parts);
double g_r(double *x, int N);

// Termodinamica
double potential_energy(
  Particles parts, LookUpTable LUT, double L, int exact);
double temperature(Particles parts);
double kinetic(Particles parts);
double pressure(
  double rho, double T, double L, LookUpTable LUT, Particles parts);

// Escribir
int write_log(
  int timestep, double time, char* filename, double rho, double L,
  LookUpTable LUT, Particles parts, int exact);
