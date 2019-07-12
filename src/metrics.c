#ifndef METRICS_H
#define METRICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metrics.h"
#include "force.h"
#include "objects.h"

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


double get_v_from_table(double r2, LookUpTable LUT) {
  /* Interpolate potential value from table.*/
  int index = (int)(LUT.length * r2 / LUT.r_c / LUT.r_c);
  if (index >= LUT.length - 1) return 0;
  else {
    double r1, r0, v1, v0;
    r1 = *(LUT.r2 + index + 1);
    r0 = *(LUT.r2 + index + 0);
    v1 = *(LUT.v + index + 1);
    v0 = *(LUT.v + index + 0);
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


double potential_energy(Particles parts, LookUpTable LUT, double L, int exact) {
  int particle_i, particle_j;
  double r2, potential;

  potential = 0;

  for (int i = 0; i < parts.N - 1; i++) {
    for (int j = i + 1; j < parts.N; j++) {
      particle_i = 3 * i;
      particle_j = 3 * j;
      r2 = r_squared(parts.x + particle_i, parts.x + particle_j, L);
      if (exact) potential += get_v_from_r2(r2, LUT.r_c);
      else potential += get_v_from_table(r2, LUT);
    }
  }
  return potential / parts.N;
}


double kinetic(Particles parts) {
  double K = 0;
  for (int i = 0; i < 3 * parts.N; i++) K += *(parts.v + i) * *(parts.v + i) / 2;
  return K / parts.N;
}


double temperature(Particles parts) {
  return 2 * kinetic(parts) / 3;
}


double pressure(double rho, double T, double L, LookUpTable LUT, Particles parts, int exact) {
  double p, r2, force, distance;
  int amount = 0;

  p = 0;
  for (int i = 0; i < parts.N - 1; i++) {
    for (int j = i + 1; j < parts.N; j++) {
      r2 = r_squared(parts.x + 3 * i, parts.x + 3 * j, L);
      if (exact) force = get_force_from_r2(r2, LUT.r_c);
      else force = get_force_from_table(r2, LUT);
      for (int dir = 0; dir < 3; dir++) {
        distance = min_diff(*(parts.x + 3 * i + dir), *(parts.x + 3 * j + dir), L);
        p += force * distance * distance;
      }
      amount++;
    }
  }
  p = p / amount / 3 / L / L / L + rho * T;
  return p;
}


double verlet_coeff(Particles parts, double L) {
  /* Get Verlet Coefficient for the disposition of particles.

  The Verlet coefficient is a meassure of the disorder of the particles.
  If they are in a cubic arrangement, this will be 1. As it get messier,
  it starts going to 0.
  */
  double lambda = 0;

  double m = L / cbrt(parts.N / 1.0);  // Characteristic distance between particles.

  for (int dir = 0; dir < 3; dir++) {
    double aux = 0;
    for (int i = 0; i < parts.N; i++) aux += cos(2 * PI / m * (parts.x[3 * i + dir] - m / 2));
    lambda += aux / parts.N;
  }
  return lambda / 3;
}


double h_boltzmann(Particles parts) {
  /* Get H from momentum distribution.

  This is a way to meassure of the entropy of the system. (kind of)
  Amount of bins is decided with Sturges Rule.
    * bins = 1 + 3. 322 logN
  */
  // int bins = (int)(1 + 3.322 * log((double)parts.N));
  int bins = 50;
  int i;
  double step;

  double* f = (double*)malloc(2 * bins * sizeof(double));
  double* v_dir = (double*)malloc(parts.N * sizeof(double));

  double H = 0;
  double Hdir = 0;

  for (int dir = 0; dir < 3; dir++) {

    for (i = 0; i < parts.N; i++) v_dir[i] = parts.v[3 * i + dir];

    double min = 1000;
    double max = -1000;
    double curr;

    for (i = 0; i < parts.N; i++) {
      curr = v_dir[i];
      if (curr < min) min = curr;
      else if (curr > max) max = curr;
    }
    step = histogram(f, v_dir, parts.N, min, max, bins);
    for (i = 0; i < bins; i++) if (*(f + i)) H += *(f + i) * log(*(f + i)) * step;

    H += Hdir / 3;
  }
  free(f);
  free(v_dir);
  return H;
}


int write_log(
    int timestep, double time, char* filename, double rho, double L,
    LookUpTable LUT, Particles parts, int exact) {
  /* Write log line with relevant metrics.*/
  FILE *fp;

  double K = kinetic(parts);
  double T = temperature(parts);
  double V = potential_energy(parts, LUT, L, exact);
  double E = V + K;
  double P = pressure(rho, T, L, LUT, parts, exact);
  double verlet = verlet_coeff(parts, L);
  double H = h_boltzmann(parts);

  if (timestep) fp = fopen(filename, "a");
  else {
    fp = fopen(filename, "w");
    fprintf(fp, "timestep,time,L,N,rho,V,K,E,T,P,verlet,H\n");
  };
  fprintf(
    fp, "%d,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
    timestep, time, L, parts.N, rho, V, K, E, T, P, verlet, H);
  fclose(fp);
  return 0;
}


#endif
