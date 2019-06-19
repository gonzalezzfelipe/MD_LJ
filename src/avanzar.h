#include "general.h"
#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>


int update_postitions(float* x, float* v, float* f, int N, float dt);
int update_velocities(float* v, float* f, int N, float dt);
int update_forces();
int timestep(float* x, float* v, float* f, int N, float dt);
