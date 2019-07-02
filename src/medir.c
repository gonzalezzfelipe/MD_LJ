#ifndef MEDIR_H
#define MEDIR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "medir.h"
#include "interaccion.h"

#define PI 3.141592


double histogram(double *y, double *x, int n, double a, double b, int m) {
  /* Generate a histogram from data.

  Get data in vector 'x' of length n and bins it in histogram 'y' of
  length 'm', in between the limits [a, b].

     y[0]...y[m-1] ==> Normalized histogram values.
     y[m]...y[2m-1] ==> mid value of each bin
  */

  int i, j;
  double h, hh, s;

  s = 1.0 / (double)n;
  h = (b - a) / (double)m;
  hh = h / 2.0;

  for (i = 0; i < m; i++) {
    *(y + i) = 0.0;
    *(y + m + i) = (double)i * h + a + hh;
  }

  for(i = 0; i < n; i++) {
    j = (int)floor((*(x + i) - a) / h);
    if (j < 0) j = 0;
    else if (j > m) j = m - 1;
    y[j] += s;
  }

  return s;
}


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


double kinetic(double *v, int N) {
  double temp = 0;
  for (int i = 0; i < 3 * N; i++) temp += *(v + i) * *(v + i) / 2;
  return temp / N;
}


double temperature(double *v, int N) {
  return kinetic(v, N) / 3;
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
  /* Get Verlet Coefficient for the disposition of particles.

  The Verlet coefficient is a meassure of the disorder of the particles.
  If they are in a cubic arrangement, this will be 1. As it get messier,
  it starts going to 0.
  */
  double lambda = 0;

  double m = L / cbrt(N / 1.0);  // Characteristic distance between particles.

  for (int dir = 0; dir < 3; dir++) {
    double aux = 0;
    for (int i = 0; i < N; i++) aux += cos(2 * PI / m * (*(x + 3 * i + dir) - m / 2));
    lambda += aux / N;
  }
  return lambda / 3;
}


double h_boltzmann(double *v, int N) {
  /* Get H from momentum distribution.

  This is a way to meassure of the entropy of the system. (kind of)
  Amount of bins is decided with Sturges Rule.

    * bins = 1 + 3. 322 logN
  */
  int bins = (int)(1 + 3.322 * log((double)N));
  int dir, i;
  double aux;

  double* p = (double*)malloc(N * sizeof(double));
  double* f = (double*)malloc(2 * bins * sizeof(double));

  double T = temperature(v, N);
  double coeff = 1 / pow(4 * PI * T, 3 / 2);

  // printf("p\n");
  for (i = 0; i < N; i++) {
    aux = 0;
    for (dir = 0; dir < 3; dir++) aux += *(v + 3 * i + dir) * *(v + 3 * i + dir);
    *(p + i) = coeff * exp(- aux / 2 / T);
    // printf("%lf\n", *(p + i));
  }

  double min = 1000;
  double max = -1;

  for (i = 0; i < N; i++) {
    if (*(p + i) < min) min = *(p + i);
    else if (*(p + i) > max) max = *(p + i);
  }

  double step = histogram(f, p, N, min, max, bins);

  double H = 0;
  for (i = 0; i < bins; i++) H += *(f + i) * log(*(f + i)) * step;
  return H;
}


int write_log(
    int timestep,
    double time,
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

  double K = kinetic(v, N);
  double T = temperature(v, N);
  double V = potential_energy(x, v, table_r2, table_v, N, L, r_c, length);
  double E = V + K;
  double P = pressure(rho, T, L, table_f, table_r2, r_c, length, x, N);
  double verlet = verlet_coeff(x, L, N);
  double H = h_boltzmann(v, N);

  if (timestep) fp = fopen(filename, "a");
  else {
    fp = fopen(filename, "w");
    fprintf(fp, "timestep,time,V,K,E,T,P,rho,verlet,H\n");
  };
  fprintf(
    fp, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
    timestep, time, V, K, E, T, P, rho, verlet, H);
  fclose(fp);
  return 0;
}


#endif
