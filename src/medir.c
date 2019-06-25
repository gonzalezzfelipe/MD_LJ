#ifndef MEDIR_H
#define MEDIR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "medir.h"
#include "interaccion.h"


float get_v_from_table(float r2, float r_c, float *table_v, float *table_r2, int length) {
  /* Interpolate potential value from table.*/
  int index = (int)(length * r2 / r_c / r_c);
  if (index >= length - 1) return 0;
  else {
    float r1, r0, f1, f0;
    r1 = *(table_r2 + index + 1);
    r0 = *(table_r2 + index + 0);
    v1 = *(table_v + index + 1);
    v0 = *(table_v + index + 0);
    return (v1 - v0) * (r2 - r0) / (r1 - r0) + v0;
  }
}

float potential_energy(float *x, float *v, float *table_r2, float *table_v, int N, int L, float r_c){
  int particle_i, particle_j;
  float r2, potential;

  for (int i = 0; i < N - 1; i++) {
    for (int j = i; j < N; j++) {
      particle_i = 3 * i;
      particle_j = 3 * j;
      r2 = r_squared(x + particle_i, x + particle_j, L);
      potential += get_v_from_table(r2, r_c, table_v, table_r2, length);
    }
  }
  return potential;
}





#endif
