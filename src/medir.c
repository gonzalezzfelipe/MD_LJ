#ifndef MEDIR_H
#define MEDIR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "medir.h"
#include "interaccion.h"

#define PI 3.141592


double get_v_from_table(double r2, double r_c, double *table_v, double *table_r2, int length) {
  /* Interpolate potential value from table.*/
  int index = (int)(length * r2 / r_c / r_c);
  if (index >= length - 1) return 0;
  else {
    double r1, r0, v1, v0;
    r1 = *(table_r2 + index + 1);
    r0 = *(table_r2 + index + 0);
    v1 = *(table_v + index + 1);
    v0 = *(table_v + index + 0);
    return (v1 - v0) * (r2 - r0) / (r1 - r0) + v0;
  }
}


double get_v_from_r2(double r2, double r_c) {
  /* Interpolate potential value from table.*/
  if (r2 > r_c * r_c) return 0;
  else {
    double r6 = r2 * r2 * r2;
    double r_c6 = r_c * r_c * r_c * r_c * r_c * r_c;
    return 4.0 * (1.0 / r6 / r6  - 1.0 / r6) - 4.0 * (1.0 / r_c6 / r_c6 - 1.0 / r_c6);
  }
}


double potential_energy(double *x, double *v, double *table_r2, double *table_v, int N, double L, double r_c, int length){
  int particle_i, particle_j;
  double r2, potential;

  potential = 0;

  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      particle_i = 3 * i;
      particle_j = 3 * j;
      r2 = r_squared(x + particle_i, x + particle_j, L);
      // potential += get_v_from_table(r2, r_c, table_v, table_r2, length);
      potential += get_v_from_r2(r2, r_c);
    }
  }
  return potential / N;
}


double temperature(double *v, int N) {
  double temp = 0;
  for (int i = 0; i < 3 * N; i++) temp += *(v + i) * *(v + i) / 2;
  return temp / N;
}


double pressure(double rho, double T, double L, double *table_f, double *table_r2, double r_c, int length, double *x, int N) {
  double p, r2, force, distance;
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


double verlet_coeff(double *x, double L, int N) {
  double lambda = 0;
  double aux = 0;

  for (int dir = 0; dir < 3; dir++) {
    for (int i = 0; i < N; i++) aux += cos(2 * PI / L * (*(x + 3 * i + dir) - L / 2));
    lambda += aux / N;
  }
  return lambda;
}


int write_log(
    int timestep,
    char* filename,
    double rho,
    int N,
    double L,
    double r_c,
    double *table_r2,
    double *table_v,
    double *table_f,
    int length,
    double *x,
    double *v,
    double *f
  ) {
  /* Write log line with relevant metrics.*/
  FILE *fp;

  double T = temperature(v, N);
  double V = potential_energy(x, v, table_r2, table_v, N, L, r_c, length);
  double E = V + T;
  double P = pressure(rho, T, L, table_f, table_r2, r_c, length, x, N);
  double verlet = verlet_coeff(x, L, N);

  if (timestep) fp = fopen(filename, "a");
  else {
    fp = fopen(filename, "w");
    fprintf(fp, "timestep,V,T,E,P,rho,verlet\n");
  };
  fprintf(fp, "%d,%lf,%lf,%lf,%lf,%lf,%lf\n", timestep, V, T, E, P, rho, verlet);
  fclose(fp);
  return 0;
}


#endif
