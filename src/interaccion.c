#ifndef INTERACCION_H
#define INTERACCION_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interaccion.h"


float min_diff(float x, float y, int periodicity) {
  /* Minimum distance between two points on a 1D periodical axis. */
  float delta;
  delta = x - y;
  if (fabsf(delta) > periodicity * 1.0 / 2)  {
    if (delta > 0) delta = periodicity - delta;
    else delta = periodicity + delta;
  }
  return delta;
}


float r_squared(float* x_1, float* x_2, int size) {
  /* Calculate the min squared distance between two particles.

  Given that the system is periodic, the distance between two points is
  not trivially the euclidean distance between them. Instead, I want the
  minimum distance between the two points and all the images.

  Returns
  -------
  float: r ** 2, the distance between the particles, squared.
  */
  int i;
  float delta, distance;

  distance = 0;
  for (i = 0; i < 3; i++) {
    delta = min_diff(*(x_1 + i), *(x_2 + i), size);
    if (delta > size * 1.0 / 2) delta = size - delta;
    distance += delta * delta;
  }
  return distance;
}


float get_force_from_table(float r2, float r_c, float *table_f, float *table_r2, int length) {
  /* Interpolate force value from table.*/
  int index = (int)(length * r2 / r_c / r_c);
  if (index >= length - 1) return 0;
  else {
    float r1, r0, f1, f0;
    r1 = *(table_r2 + index + 1);
    r0 = *(table_r2 + index + 0);
    f1 = *(table_f + index + 1);
    f0 = *(table_f + index + 0);
    return (f1 - f0) * (r2 - r0) / (r1 - r0) + f0;
  }
}


int update_forces(float *f, float *x, int N, int L, float r_c, float *table_f, float *table_r2, int length) {
  /* Update forces vector using table. */
  float distance;
  float force;

  for (int i = 0; i < N - 1; i++) {
    *(f + 3 * i) = 0;
    *(f + 3 * i + 1) = 0;
    *(f + 3 * i + 2) = 0;
    for (int j = i + 1; j < N; j++) {
      distance = r_squared(x + 3 * i, x + 3 * j, L);
      if (distance < r_c * r_c) {
        force = get_force_from_table(distance, r_c, table_f, table_r2, length);
        for (int dir = 0; dir < 3; dir++) {
          *(f + 3 * i + dir) += force * (*(x + 3 * i + dir) - *(x + 3 * j + dir));
          *(f + 3 * j + dir) -= force * (*(x + 3 * i + dir) - *(x + 3 * j + dir));
        }
      }
    }
  }
  return 0;
}


#endif
