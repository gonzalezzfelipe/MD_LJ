#include "interaccion.h"
#include "init.h"
#include "avanzar.h"
#include "objetos.h"
#include "visualizacion.h"
#include "medir.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define DEF_DT 0.001
#define DEF_FRAMES 100
#define DEF_FRAMES_STEP 1
#define DEF_LAMMPSTRJ_FILENAME "data/test.lammpstrj"
#define DEF_LOG_FILENAME "data/test.log"
#define DEF_TERMALIZATION 2500
#define R_C 2.5
#define TABLE_LENGTH 10000
#define RESCALING_RELATIVE_ERROR 1.01
#define INITIAL_TEMPERATURE 4.0


int main(int argc, char *argv[]){
  /* Run MD simulation.

  Parameters
  ----------
  double N:
    Amount of particles.
  double rho:
    Desired density.
  double T:
    Desired temperature.
  double dt, optional, default 0.001:
    Temporal step.
  int frames, optional, default 100:
    Amount of frames to save on lammpstrj file.
  int frames_step, optional, default 1:
    Amount of frames to leave between saved frames.
  int termalization, optional, default 2500:
    Amount of timesteps for termalization.
  char lamppstrj_filename, optional, default: data/test.lammpstrj
    Filename where to write results.
  char log_filename, optional, default: data/test.lammpstrj
    Filename where to write results.
  int seed, optional
    Seed to feed random number generator. Default is current timestamp.
  */
  int frames, frames_step, N, termalization;
  int* seed;
  double T, rho, dt;
  char lamppstrj_filename[255];
  char log_filename[255];

  seed = (int*)malloc(sizeof(int));

  sscanf(argv[1], "%d", &N);
  sscanf(argv[2], "%lf", &rho);
  sscanf(argv[3], "%lf", &T);
  if (argc >= 5) sscanf(argv[4], "%lf", &dt);
  else dt = 0.001;
  if (argc >= 6) sscanf(argv[5], "%d", &frames);
  else frames = DEF_FRAMES;
  if (argc >= 7) sscanf(argv[6], "%d", &frames_step);
  else frames_step = DEF_FRAMES_STEP;
  if (argc >= 8) sscanf(argv[7], "%d", &termalization);
  else termalization = DEF_TERMALIZATION;
  if (argc >= 9) strcpy(lamppstrj_filename, argv[8]);
  else strcpy(lamppstrj_filename, DEF_LAMMPSTRJ_FILENAME);
  if (argc >= 10) strcpy(log_filename, argv[9]);
  else strcpy(log_filename, DEF_LOG_FILENAME);
  if (argc >= 11) sscanf(argv[10], "%d", &*seed);
  else *seed = rand();

  int L;

  // Initialize
  L = cbrt(N / 1.0 / rho);

  struct Particles particles;

  particles.x = (double*)malloc(3 * N * sizeof(double));
  particles.v = (double*)malloc(3 * N * sizeof(double));
  particles.f = (double*)malloc(3 * N * sizeof(double));
  particles.N = N;

  struct LookUpTable LUT;

  LUT.r = (double*)malloc(TABLE_LENGTH * sizeof(double));
  LUT.r2 = (double*)malloc(TABLE_LENGTH * sizeof(double));
  LUT.f = (double*)malloc(TABLE_LENGTH * sizeof(double));
  LUT.v = (double*)malloc(TABLE_LENGTH * sizeof(double));
  LUT.length = TABLE_LENGTH;
  LUT.r_c = R_C;

  rho = initial_positions(L, particles);
  initial_velocities(particles, INITIAL_TEMPERATURE);

  fill_forces_table(LUT);
  update_forces(particles, L, LUT);

  // Termalize
  if (termalization) rescaling(T, RESCALING_RELATIVE_ERROR, termalization, particles, dt, L, LUT);

  // Evolve
  double time = 0;
  for (int frame = 0; frame < frames; frame++) {
    time = frame * frames_step * dt;
    save_lammpstrj(lamppstrj_filename, particles.x, particles.v, N, L, frame);
    write_log(frame, time, log_filename, rho, L, LUT, particles);
    for (int i = 0; i < frames_step; i++) timestep(particles, dt, L, LUT);
  }

  free(particles.x);
  free(particles.v);
  free(particles.f);
  free(LUT.r);
  free(LUT.r2);
  free(LUT.v);
  free(LUT.f);
  return 0;
}
