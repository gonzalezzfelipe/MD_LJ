#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>


int update_postitions(float* x, float* v, float* f, int N, float dt, int L);
int update_velocities(float* v, float* f, int N, float dt);
int timestep(float* x, float* v, float* f, int N, float dt, int L, float r_c, float *table_f, float *table_r2, int length);
