#ifndef AVANZAR_H
#define AVANZAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "avanzar.h"
#include "interaccion.h"
#include "objetos.h"


int update_postitions(Particles parts, double dt, double L) {
  /* Calculate new positions for the particles.

  Uses Verlet Algorithm:

    *) x(t + dt) = x(t) + v(t) * dt + 1 / 2 * f / m * dt ** 2

  Assumes m = 1.
  */
  for (int i = 0; i < 3 * parts.N; i++) {
    parts.x[i] += parts.v[i] * dt + parts.f[i] * dt * dt / 2;
    if (parts.x[i] < 0) parts.x[i] += L;
    else if (parts.x[i] > L) parts.x[i] -= L;
  }
  return 0;
}


int update_velocities(Particles parts, double dt) {
  /* Calculate intermediate step for new velocities for the particles.

  Uses Verlet Algorithm:

    *) v(t + h) = v(t) + dt / 2 * (f(t + h) + f(t))

  And it breaks it down into 2 v += dt / 2m * f, this way the forces in
  both times dont have to be known simultaneously.

  Assumes m = 1.
  */
  for (int i = 0; i < 3 * parts.N; i++) parts.v[i] += parts.f[i] * 0.5 * dt;
  return 0;
}


int timestep(Particles parts, double dt, double L, LookUpTable LUT) {
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
  update_postitions(parts, dt, L);
  update_velocities(parts, dt);
  update_forces(parts, L, LUT);
  update_velocities(parts, dt);
  return 0;
}


#endif
