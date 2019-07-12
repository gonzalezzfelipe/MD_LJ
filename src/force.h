#include <stdio.h>
#include <stdlib.h>
#include "objects.h"


double min_diff(double x, double y, double L);
double r_squared(double* x_1, double* x_2, double L);
double potential(double r_squared, double sigma);
double get_force_from_table(double r2, LookUpTable LUT);
double get_force_from_r2(double r2, double r_c);
int update_forces(Particles parts, double L, LookUpTable LUT, int exact);
