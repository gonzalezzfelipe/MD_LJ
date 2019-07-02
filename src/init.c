#ifndef INIT_H
#define INIT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"
#include "medir.h"
#include "avanzar.h"

#define PI 3.141592


double sample_boltzmann(double T) {
  /* Get random sample from Maxwell - Boltzmann distribution.

  To sample this value, 10 different random samples between [-1, 1] and
  sum them up, as this is equal to a sample from a random distribution
  of mean 0 and variance 10 / 3.

  The velocities are defined following a Maxwell - Boltzmann distribution.

    *) dN / N = (2 * pi * kb * T / m) ^ (-1 / 2) exp(- m v ** 2 / 2 kb T)

  m and kb are assumed to be 1.

  Parameters
  ----------
  double T:
    T of the sample.

  Returns
  -------
  double: Random number sampled from Boltzmann distribution.
  */
  int i;
  double sample;

  sample = 0.0;
  for(i = 0; i < 10; i++) sample += 2 * rand() / 1.0 / RAND_MAX - 1;
  sample = 3 * T * sample / 10;
  return sample;
}


int amount_of_particles(double rho, double L) {
  /* Define amount of particles in arangement.

  Given a density rho and a volume L ** 3, the amount of particles are
  defined as the amount of particles that fit on that volume.

  Parameters
  ----------
  double rho:
    Desired density of the arangement.
  double L:
    Side of the cubic volume where to place the particles.

  Returns
  -------
  int:
    Amount of particles.
  */
  int n_x;
  n_x = (int)(cbrt(rho) * L);
  return n_x * n_x * n_x;
}


double initial_positions(double L, double* x, int n) {
  /* Define initial positions for all particles.

  The particles will be set up on a simple qubic arangement.

  Parameters
  ----------
  double L:
    Side of the cubic volume where to place the particles.
  double* x:
    Pointer to the vector that has the positions of every particle as:
      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  int n:
    Amount of particles.
  */
  int i, j, k, n_x, particle;
  double step;

  n_x = (int)cbrt(n);
  step = L / n_x;

  for (i = 0; i < n_x; i++) {
    for (j = 0; j < n_x; j++) {
      for (k = 0; k < n_x; k++) {
        particle = n_x * n_x * i + n_x * j + k;
        *(x + 3 * particle + 0) = step / 2 + step * i;  // x
        *(x + 3 * particle + 1) = step / 2 + step * j;  // y
        *(x + 3 * particle + 2) = step / 2 + step * k;  // z
      }
    }
  }
  return n / 1.0 / L / L / L;
}


int initial_velocities(double* v, int n, double T) {
  /* Fill initial velocities vector.

  The velocities are defined following a Maxwell - Boltzmann distribution.

    *) dN / N = (2 * pi * kb * T / m) ^ (-1 / 2) exp(- m v ** 2 / 2 kb T)

  Given a T, a random sample is taken from this distribution for
  each particle direction. Finally the mean velocity is substracted in each
  direction as for the sample not to move uniformly in a direction.

  Parameters
  ----------
  double* v:
    Pointer to the vector that has the velocities of every particle as:
      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  int n:
    Amount of particles.
  double T:
    T of the sample.
  */
  int i, direction;

  double mean;

  for (i = 0; i < 3 * n; i++) *(v + i) = sample_boltzmann(T);
  for (direction = 0; direction < 3; direction++) {
    mean = 0;
    for (i = 0; i < n; i++) mean += *(v + 3 * i + direction) / n;
    for (i = 0; i < n; i++) *(v + 3 * i + direction) -= mean;
  }
  return 0;
}


int fill_forces_table(
    double *table_r,
    double *table_r2,
    double *table_f,
    double *table_v,
    double r_c,
    int length) {
  /* Create table including r, r squared, F and V.

  The idea is to have a table, indexed by r squared as not to have to
  calculate the force at each step, and replace it with a lookup at this
  table plus an interpolation.

  NOTE: Indexed forces correspond to F / r, as to be able to do:
    * Fx = x * table_f,
    * Fy = y * table_f,
    * Fz = z * table_f
  */
 double step = r_c * r_c / length;
 double current = 0.0;
 double current6;

 double v_c = 4.0 / pow(r_c, 12) - 4.0 / pow(r_c, 6);

 for (int i = 1; i < length; i++) {
   current += step;

   *(table_r + i) = sqrt(current);
   *(table_r2 + i) = current;

   current6 = current * current * current;
   *(table_f + i) = 48.0 / current6 / current6 / current - 24.0 / current6 / current;
   *(table_v + i) = 4.0 / current6 / current6 - 4.0 / current6 - v_c;
 }
 return 0;
}

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
    int length) {
  /* Apply reescaling until desired temperature is achieved. */
  int i;
  float coeff;
  double actual = temperature(v, N);
  double error = fabs(actual - T) / T;

  while (error > relative_error) {
    for (i = 0; i < termalization; i++) timestep(x, v, f, N, dt, L, r_c, table_f, table_r2, length);
    actual = temperature(v, N);
    coeff = sqrt(T / actual);
    for (i = 0; i < 3 * N; i++) *(v + i) = coeff * *(v + i);
    printf("Desired Temperature: %f\n", T);
    printf("Actual Temperature: %f\n", actual);
  }
  return 0;
}


#endif
