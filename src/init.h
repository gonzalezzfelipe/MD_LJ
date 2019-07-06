#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "objetos.h"

#define PI 3.141592


double sample_boltzmann(double temperature);
int amount_of_particles(double rho, double L);
double initial_positions(double L, Particles parts);
int initial_velocities(Particles parts, double temperature);
int fill_lut(LookUpTable LUT, float r_c, int length);
int init_particles(Particles *particles, int N, float initial_t, float L, LookUpTable LUT);

int rescaling(
  double T, double relative_error, int termalization, Particles parts,
  double dt, double L, LookUpTable LUT, int verbose);
