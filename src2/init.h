#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592


float sample_boltzmann(float temperature);
int amount_of_particles(float rho, int size);
float initial_positions(int size, float* x, int n);
int initial_velocities(float* v, int n, float temperature);


float sample_boltzmann(float T) {
  /* Get random sample from Maxwell - Boltzmann distribution.

  To sample this value, 10 different random samples between [-1, 1] and
  sum them up, as this is equal to a sample from a random distribution
  of mean 0 and variance 10 / 3.

  The velocities are defined following a Maxwell - Boltzmann distribution.

    *) dN / N = (2 * pi * kb * T / m) ^ (-1 / 2) exp(- m v ** 2 / 2 kb T)

  m and kb are assumed to be 1.

  Parameters
  ----------
  float T:
    T of the sample.

  Returns
  -------
  float: Random number sampled from Boltzmann distribution.
  */
  int i;
  float sample;

  sample = 0.0;
  for(i = 0; i < 10; i++) sample += 2 * rand() / 1.0 / RAND_MAX - 1;
  sample = 3 * T * sample / 10;
  return sample;
}


int amount_of_particles(float rho, int L) {
  /* Define amount of particles in arangement.

  Given a density rho and a volume L ** 3, the amount of particles are
  defined as the amount of particles that fit on that volume.

  Parameters
  ----------
  float rho:
    Desired density of the arangement.
  int L:
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


float initial_positions(int L, float* x, int n) {
  /* Define initial positions for all particles.

  The particles will be set up on a simple qubic arangement.

  Parameters
  ----------
  int L:
    Side of the cubic volume where to place the particles.
  float* x:
    Pointer to the vector that has the positions of every particle as:
      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  int n:
    Amount of particles.
  */
  int i, j, k, n_x, particle;
  float step;

  n_x = (int)cbrt(n);
  step = L / 1.0 / n_x;

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


int initial_velocities(float* v, int n, float T) {
  /* Fill initial velocities vector.

  The velocities are defined following a Maxwell - Boltzmann distribution.

    *) dN / N = (2 * pi * kb * T / m) ^ (-1 / 2) exp(- m v ** 2 / 2 kb T)

  Given a T, a random sample is taken from this distribution for
  each particle direction. Finally the mean velocity is substracted in each
  direction as for the sample not to move uniformly in a direction.

  Parameters
  ----------
  float* v:
    Pointer to the vector that has the velocities of every particle as:
      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  int n:
    Amount of particles.
  float T:
    T of the sample.
  */
  int i, direction;

  float mean;

  for (i = 0; i < 3 * n; i++) *(v + i) = sample_boltzmann(T);
  for (direction = 0; direction < 3; direction++) {
    mean = 0;
    for (i = 0; i < n; i++) mean += *(v + 3 * i + direction) / n;
    for (i = 0; i < n; i++) *(v + 3 * i + direction) -= mean;
  }
  return 0;
}
