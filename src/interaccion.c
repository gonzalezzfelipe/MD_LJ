#ifndef INTERACCION_H
#define INTERACCION_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interaccion.h"
#include "objetos.h"


double min_diff(double x, double y, double L) {
  /* Minimum distance between two points on a 1D periodical axis. */
  double delta;
  delta = x - y;
  if (delta > L / 2 || delta < - L / 2) {
    if (delta > 0) delta = delta - L;
    else delta = L + delta;
  }
  return delta;
}


double r_squared(double* x_1, double* x_2, double L) {
  /* Calculate the min squared distance between two particles.

  Given that the system is periodic, the distance between two points is
  not trivially the euclidean distance between them. Instead, I want the
  minimum distance between the two points and all the images.

  Returns
  -------
  double: r ** 2, the distance between the particles, squared.
  */
  int i;
  double delta, distance;

  distance = 0;
  for (i = 0; i < 3; i++) {
    delta = min_diff(*(x_1 + i), *(x_2 + i), L);
    distance += delta * delta;
  }
  return distance;
}


double get_force_from_table(double r2, LookUpTable LUT) {
  /* Interpolate force value from table.*/
  int index = (int)(LUT.length * r2 / LUT.r_c / LUT.r_c);
  if (index >= LUT.length - 1) return 0;
  else {
    double r1, r0, f1, f0;
    r1 = *(LUT.r2 + index + 1);
    r0 = *(LUT.r2 + index + 0);
    f1 = *(LUT.f + index + 1);
    f0 = *(LUT.f + index + 0);
    return (f1 - f0) * (r2 - r0) / (r1 - r0) + f0;
  }
}


double get_force_from_r2(double r2, double r_c) {
  /* Get Lennard Jones force from r squared. */
  if (r2 > r_c * r_c) return 0;
  else {
    double r6 = r2 * r2 * r2;
    return 48.0 / r6 / r6 / r2 - 24.0 / r6 / r2;
  }
}



int update_forces(Particles parts, double L, LookUpTable LUT) {
  /* Update forces vector using table. */
  double distance, r2;
  double force;

  int i, j, dir;

  for (i = 0; i < 3 * parts.N; i++) parts.f[i] = 0;

  for (i = 0; i < parts.N - 1; i++) {
    for (j = i + 1; j < parts.N; j++) {
      r2 = r_squared(parts.x + 3 * i, parts.x + 3 * j, L);
      if (r2 < LUT.r_c * LUT.r_c) {
        force = get_force_from_table(r2, LUT);
        // force = get_force_from_r2(r2, r_c);
        for (dir = 0; dir < 3; dir++) {
          distance = min_diff(parts.x[3 * i + dir], parts.x[3 * j + dir], L);
          parts.f[3 * i + dir] += force * distance;
          parts.f[3 * j + dir] -= force * distance;
        }
      }
    }
  }
  return 0;
}


#endif
