#ifndef AVANZAR_H
#define AVANZAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "avanzar.h"
#include "interaccion.h"


int update_postitions(double* x, double* v, double* f, int N, double dt, double L) {
  /* Calculate new positions for the particles.

  Uses Verlet Algorithm:

    *) x(t + dt) = x(t) + v(t) * dt + 1 / 2 * f / m * dt ** 2

  Assumes m = 1.

  Parameters
  ----------
  double* x:
    Pointer to the vector that has the positions of every particle as:

      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  double* v:
    Pointer to the vector that has the velocities of every particle as:

      v = [v_x_0, v_y_0, v_z_0, v_x_1, v_y_1, v_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  double* f:
    Pointer to the vector that has the forces of every particle as:

      f = [f_x_0, f_y_0, f_z_0, f_x_1, f_y_1, f_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  int N:
    Amount of particles. The length of the other vectors should be 3 * n
  double dt:
    Temporal step.
  double L
    Size of the box.
  */
  for (int i = 0; i < 3 * N; i++) {
    *(x + i) += *(v + i) * dt + *(f + i) * dt * dt / 2;
    if (*(x + i) < 0) *(x + i) += L;
    else if (*(x + i) > L) *(x + i) -= L;
  }
  return 0;
}


int update_velocities(double* v, double* f, int N, double dt) {
  /* Calculate intermediate step for new velocities for the particles.

  Uses Verlet Algorithm:

    *) v(t + h) = v(t) + dt / 2 * (f(t + h) + f(t))

  And it breaks it down into 2 v += dt / 2m * f, this way the forces in
  both times dont have to be known simultaneously.

  Assumes m = 1.

  Parameters
  ----------
  double* x:
    Pointer to the vector that has the positions of every particle as:

      x = [x_0, y_0, z_0, x_1, y_1, z_1, ...]
          --------------  -------------
            Particle 0     Particle 1
  double* v:
    Pointer to the vector that has the velocities of every particle as:

      v = [v_x_0, v_y_0, v_z_0, v_x_1, v_y_1, v_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  double* f:
    Pointer to the vector that has the forces of every particle as:

      f = [f_x_0, f_y_0, f_z_0, f_x_1, f_y_1, f_z_1, ...]
          --------------------  -------------------
                Particle 0          Particle 1
  int N:
    Amount of particles. The length of the other vectors should be 3 * n
  double dt:
    Temporal step.
  */
  for (int i = 0; i < 3 * N; i++) *(v + i) = *(v + i) + *(f + i) * 0.5 * dt;
  return 0;
}


int timestep(double* x, double* v, double* f, int N, double dt, double L, double r_c, double *table_f, double *table_r2, int length) {
  /* Make a step in time.

  This includes updating positions, velocities and forces using the
  Verlet algorithm:

    * x(t + h) = x(t) + v(t) * dt + dt ** 2 / 2m f(t)
    * v(t + h) = v(t) + dt / 2m * (f(t) + f(t + h))

  The step in velocities depends on the forces after and before the
  step, this can be calculated because forces depend solely on x.

  To implement this, first new positions are calculated, then:

    v = v + dt / 2m f
    update_forces()
    v = v + dt / 2m f

  This way there's no need to know the forces at both times simultaneously.
  */
  update_postitions(x, v, f, N, dt, L);
  update_velocities(v, f, N, dt);
  update_forces(f, x, N, L, r_c, table_f, table_r2, length);
  update_velocities(v, f, N, dt);
  return 0;
}


#endif
