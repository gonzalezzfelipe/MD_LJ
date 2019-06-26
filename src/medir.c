#ifndef MEDIR_H
#define MEDIR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "medir.h"
#include "interaccion.h"

#define PI 3.141592


float get_v_from_table(float r2, float r_c, float *table_v, float *table_r2, int length) {
  /* Interpolate potential value from table.*/
  int index = (int)(length * r2 / r_c / r_c);
  if (index >= length - 1) return 0;
  else {
    float r1, r0, v1, v0;
    r1 = *(table_r2 + index + 1);
    r0 = *(table_r2 + index + 0);
    v1 = *(table_v + index + 1);
    v0 = *(table_v + index + 0);
    return (v1 - v0) * (r2 - r0) / (r1 - r0) + v0;
  }
}

float potential_energy(float *x, float *v, float *table_r2, float *table_v, int N, float L, float r_c, int length){
  int particle_i, particle_j;
  float r2, potential;

  potential = 0;

  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      particle_i = 3 * i;
      particle_j = 3 * j;
      r2 = r_squared(x + particle_i, x + particle_j, L);
      potential += get_v_from_table(r2, r_c, table_v, table_r2, length);
    }
  }
  return potential / N;
}


float temperature(float *v, int N) {
  float temp = 0;
  for (int i = 0; i < 3 * N; i++) temp += *(v + i) * *(v + i);
  return temp / 3 / N;
}


float pressure(float rho, float T, float L, float *table_f, float *table_r2, float r_c, int length, float *x, int N) {
  float p, r2, force, distance;
  int amount = 0;

  p = 0;
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      r2 = r_squared(x + 3 * i, x + 3 * j, L);
      force = get_force_from_table(r2, r_c, table_f, table_r2, length);
      for (int dir = 0; dir < 3; dir++) {
        distance = min_diff(*(x + 3 * i + dir), *(x + 3 * j + dir), L);
        p += force * distance * distance;
      }
      amount++;
    }
  }
  p = p / amount / 3 / L / L / L + rho * T;
  return p;
}


float verlet_coeff(float *x, float L, int N) {
  float lambda = 0;
  float aux = 0;

  for (int dir = 0; dir < 3; dir++) {
    for (int i = 0; i < N; i++) aux += cos(2 * PI / L * (*(x + 3 * i + dir) - L / 2));
    lambda += aux / N;
  }
  return lambda;
}


int write_log(
    int timestep,
    char* filename,
    float rho,
    int N,
    float L,
    float r_c,
    float *table_r2,
    float *table_v,
    float *table_f,
    int length,
    float *x,
    float *v,
    float *f
  ) {
  /* Write log line with relevant metrics.*/
  FILE *fp;

  float T = temperature(v, N);

  if (timestep) fp = fopen(filename, "a");
  else {
    fp = fopen(filename, "w");
    fprintf(fp, "timestep,V,T,P,rho,verlet\n");
  };
  fprintf(
    fp, "%d,%f,%f,%f,%f,%f\n",
    timestep,
    potential_energy(x, v, table_r2, table_v, N, L, r_c, length),
    T,
    pressure(rho, T, L, table_f, table_r2, r_c, length, x, N),
    rho,
    verlet_coeff(x, L, N)
  );
  fclose(fp);
  return 0;
}


#endif
