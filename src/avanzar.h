#ifndef AVANZAR_H
#define AVANZAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int update_postitions(float* x, float* v, float* f, int n, float dt) {
  /* Calculate new positions for the particles.

  Uses Verlet Algorithm:

    *) x(t + dt) = x(t) + v(t) * dt + 1 / 2 * f / m * dt ** 2

  Assumes m = 1.

  Parameters
  ----------
  float* x:
    Pointer to the vector that has the positions of every particle as:

      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  float* v:
    Pointer to the vector that has the velocities of every particle as:

      v = [v_x_0, v_y_0, v_z_0, v_x_1, v_y_1, v_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  float* f:
    Pointer to the vector that has the forces of every particle as:

      f = [f_x_0, f_y_0, f_z_0, f_x_1, f_y_1, f_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  int n:
    Amount of particles. The length of the other vectors should be 3 * n
  float dt:
    Temporal step.
  */
  int particle, direction, aux;
  for (particle = 0; particle < n; particle++) {
    for (direction = 0; direction < 3; direction++) {
      aux = 3 * particle + j;
      *(x + aux) = *(x + aux) + *(v + aux) * dt + *(2 + aux) * dt * dt / 2
    }
  }
  return 0;
}


int update_velocities(float* x, float* v, float* f, int n, float dt) {
  /* Calculate new velocities for the particles.

  Uses Verlet Algorithm:

    *) v(t + h) = v(t) + dt / 2 * (f(t + h) + f(t))

  Assumes m = 1.

  Parameters
  ----------
  float* x:
    Pointer to the vector that has the positions of every particle as:

      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  float* v:
    Pointer to the vector that has the velocities of every particle as:

      v = [v_x_0, v_y_0, v_z_0, v_x_1, v_y_1, v_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  float* f:
    Pointer to the vector that has the forces of every particle as:

      f = [f_x_0, f_y_0, f_z_0, f_x_1, f_y_1, f_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  int n:
    Amount of particles. The length of the other vectors should be 3 * n
  float dt:
    Temporal step.
  */
  return 0;
}


int update_forces(float* x, float* v, float* f, int n, float dt) {
  /* Calculate forces for each particle.

  */
}

#endif
