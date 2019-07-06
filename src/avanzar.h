#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h> 


int update_postitions(Particles parts, double dt, double L);
int update_velocities(Particles parts, double dt);
int timestep(Particles parts, double dt, double L, LookUpTable LUT);
