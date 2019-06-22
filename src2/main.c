// #include "interaccion.h"
#include "init.h"
#include "avanzar.h"
#include "visualizacion.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define DEF_DT 0.001
#define DEF_FRAMES 100
#define DEF_FRAMES_STEP 0
#define DEF_FILENAME "test.lammpstrj"


int main(int argc, char *argv[]){
  /* Run MD simulation.

  Parameters
  ----------
  int L:
    Size of box.
  float rho:
    Desired density.
  float T:
    Desired temperature.
  float dt, optional, default 0.001:
    Temporal step.
  int frames, optional, default 100:
    Amount of frames to save on lammpstrj file.
  int frames_step, optional, default 0:
    Amount of frames to leave between saved frames.
  char filename, optional, default: test.lammpstrj
    Filename where to write results.
  int seed, optional
    Seed to feed random number generator. Default is current timestamp.
  */
  int L, frames, frames_step;
  int* seed;
  float T, rho, dt;
  char filename[255];

  seed = (int*)malloc(sizeof(int));

  sscanf(argv[1], "%d", &L);
  sscanf(argv[2], "%f", &rho);
  sscanf(argv[3], "%f", &T);
  if (argc >= 5) sscanf(argv[4], "%f", &dt);
  else dt = 0.001;
  if (argc >= 6) sscanf(argv[5], "%d", &frames);
  else frames = DEF_FRAMES;
  if (argc >= 7) sscanf(argv[6], "%d", &frames_step);
  else frames_step = DEF_FRAMES_STEP;
  if (argc >= 8) strcpy(filename, argv[7]);
  else strcpy(filename, DEF_FILENAME);
  if (argc >= 9) sscanf(argv[8], "%d", &*seed);
  else *seed = rand();

  int N;

  // Initialize
  N = amount_of_particles(rho, L);
  printf("Amount of particles is: %d\n", N);

  float* x = (float*)malloc(3 * N * sizeof(float));
  float* v = (float*)malloc(3 * N * sizeof(float));
  float* f = (float*)malloc(3 * N * sizeof(float));

  rho = initial_positions(L, x, N);
  printf("Actual density is: %f\n", rho);
  initial_velocities(v, N, T);

  // Termalize

  // Evolve
  for (int frame = 0; frame < frames; frame++) {
    save_lammpstrj(filename, x, v, N, L, frame);
    for (int i = 0; i < frames_step + 1; i++) timestep(x, v, f, N, dt, L);
  }

  free(x);
  free(v);
  free(f);
  return 0;
}
