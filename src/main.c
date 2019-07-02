#include "interaccion.h"
#include "init.h"
#include "avanzar.h"
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
#define DEF_LAMMPSTRJ_FILENAME "test.lammpstrj"
#define DEF_LOG_FILENAME "test.log"
#define R_C 2.5
#define TABLE_LENGTH 10000


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
  char lamppstrj_filename, optional, default: test.lammpstrj
    Filename where to write results.
  int seed, optional
    Seed to feed random number generator. Default is current timestamp.
  */
  int frames, frames_step, N;
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
  if (argc >= 8) strcpy(lamppstrj_filename, argv[7]);
  else strcpy(lamppstrj_filename, DEF_LAMMPSTRJ_FILENAME);
  if (argc >= 9) strcpy(log_filename, argv[8]);
  else strcpy(log_filename, DEF_LOG_FILENAME);
  if (argc >= 10) sscanf(argv[9], "%d", &*seed);
  else *seed = rand();

  int L;

  // Initialize
  L = cbrt(N / 1.0 / rho);

  double* x = (double*)malloc(3 * N * sizeof(double));
  double* v = (double*)malloc(3 * N * sizeof(double));
  double* f = (double*)malloc(3 * N * sizeof(double));

  double* table_r = (double*)malloc(TABLE_LENGTH * sizeof(double));
  double* table_r2 = (double*)malloc(TABLE_LENGTH * sizeof(double));
  double* table_f = (double*)malloc(TABLE_LENGTH * sizeof(double));
  double* table_v = (double*)malloc(TABLE_LENGTH * sizeof(double));

  rho = initial_positions(L, x, N);
  printf("Actual density is: %lf\n", rho);
  initial_velocities(v, N, T);

  fill_forces_table(table_r, table_r2, table_f, table_v, R_C, TABLE_LENGTH);
  update_forces(f, x, N, L, R_C, table_f, table_r2, TABLE_LENGTH);

  // Termalize

  // Evolve
  printf("frame,r,r2,v,f\n");
  for (int frame = 0; frame < frames; frame++) {
    save_lammpstrj(lamppstrj_filename, x, v, N, L, frame);
    write_log(frame, log_filename, rho, N, L, R_C, table_r2, table_v, table_f, TABLE_LENGTH, x, v, f);
    for (int i = 0; i < frames_step; i++) timestep(x, v, f, N, dt, L, R_C, table_f, table_r2, TABLE_LENGTH);
  }

  free(x);
  free(v);
  free(f);
  return 0;
}
