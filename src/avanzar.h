#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>


int update_postitions(double* x, double* v, double* f, int N, double dt, double L);
int update_velocities(double* v, double* f, int N, double dt);
int timestep(double* x, double* v, double* f, int N, double dt, double L, double r_c, double *table_f, double *table_r2, int length);
