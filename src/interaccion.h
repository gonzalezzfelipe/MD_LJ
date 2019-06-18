#ifndef INTERACCION_H
#define INTERACCION_H

#include <stdlib.h>
#include <stdio.h>

float r_squared(float* x_1, float* x_2, int size) {
  /* Calculate the min squared distance between two particles.

  Given that the system is periodic, the distance between two points is
  not trivially the euclidean distance between them. Instead, I want the
  minimum distance between the two points and all the images.

  Returns
  -------
  float: r ** 2, the distance between the particles, squared.
  */
  int i;
  float delta, distance;

  distance = 0;
  for (i = 0; i < 3; i++) {
    delta = min_diff(*(x_1 + i), *(x_2 + i), size);
    if (delta > size * 1.0 / 2) delta = size - delta;
    distance += delta * delta;
  }
  return distance;
}


float potential(float r_squared, float sigma) {
  /* Calculate Lennard Jones force for a fixed distance.

  NOTE: The function depends assumes the value provided is r_squared and
  the return value is so that if multiplied by the direction projection
  you get the force.
   */
   float f, r_sixth;

   r_sixth = r_squared * r_squared * r_squared;

   f = 24 * (sigma / (r_sixth * r_squared));
}


float min_diff(float x, float y, int periodicity) {
  /* Minimum distance between two points on a 1D periodical axis. */
  float delta;
  delta = x - y;
  if (abs(delta) > size * 1.0 / 2) delta = size - delta;
  return delta;
}


#endif
