#include "interaccion.h"
#include "init.h"
#include "avanzar.h"
#include "objetos.h"
#include "visualizacion.h"
#include "argumentos.h"
#include "medir.h"
#include <argp.h>
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
#define DEF_R_C 2.5
#define DEF_TABLE_LENGTH 10000
#define DEF_RESCALING_RELATIVE_ERROR 1.01
#define DEF_INITIAL_TEMPERATURE 4.0


int main(int argc, char *argv[]){
  // Positional Arguments
  int N;
  double rho;
  double T;

  struct arguments arguments;

  // Default options
  arguments.dt = DEF_DT;
  arguments.frames = DEF_FRAMES;
  arguments.frames_step = DEF_FRAMES_STEP;
  arguments.termalization = DEF_TERMALIZATION;
  strcpy(arguments.lamppstrj_filename, DEF_LAMMPSTRJ_FILENAME);
  strcpy(arguments.log_filename, DEF_LOG_FILENAME);
  arguments.seed = time(NULL);
  arguments.r_c = DEF_R_C;
  arguments.table_length = DEF_TABLE_LENGTH;
  arguments.rescaling_relative_error = DEF_RESCALING_RELATIVE_ERROR;
  arguments.initial_temperature = DEF_INITIAL_TEMPERATURE;
  arguments.verbose = 0;

  // Parse the arguments
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  // Parse de positional arguments.
  sscanf(arguments.args[0], "%d", &N);
  sscanf(arguments.args[1], "%lf", &rho);
  sscanf(arguments.args[2], "%lf", &T);

  // Parameters
  if (arguments.verbose) {
    printf("Parameters\n");
    printf("==========\n");
    printf("Amount of particles: %d\n", N);
    printf("Density: %lf\n", rho);
    printf("Temperature: %lf\n", T);
    printf("Temporal step: %lf\n", arguments.dt);
    printf("Amount of measurements: %d\n", arguments.frames);
    printf("Amount of temporal steps between measurements: %d\n", arguments.frames_step);
    printf("Name of LAMPPSTRJ file: %s\n", arguments.lamppstrj_filename);
    printf("Name of log file: %s\n", arguments.log_filename);
    printf("Cut radius: %lf\n", arguments.r_c);
    printf("LUT table length: %d\n", arguments.table_length);
    printf("Rescaling relative error: %lf\n", arguments.rescaling_relative_error);
    printf("Initial temperature for rescaling: %lf\n", arguments.initial_temperature);
    printf("Random seed: %d\n", arguments.seed);
  }

  // Set random seed.
  srand(arguments.seed);

  // Initialize
  double L;
  L = cbrt(N / 1.0 / rho);

  struct LookUpTable LUT;
  struct Particles particles;

  LUT.r = (double*)malloc(arguments.table_length * sizeof(double));
  LUT.r2 = (double*)malloc(arguments.table_length * sizeof(double));
  LUT.f = (double*)malloc(arguments.table_length * sizeof(double));
  LUT.v = (double*)malloc(arguments.table_length * sizeof(double));
  LUT.r_c = arguments.r_c;
  LUT.length = arguments.table_length;

  fill_lut(LUT, arguments.r_c, arguments.table_length);

  particles.x = (double*)malloc(3 * N * sizeof(double));
  particles.v = (double*)malloc(3 * N * sizeof(double));
  particles.f = (double*)malloc(3 * N * sizeof(double));
  particles.N = N;

  initial_positions(L, particles);

  // Termalize
  if (arguments.termalization) {
    initial_velocities(particles, arguments.initial_temperature);
    rescaling(
      T, arguments.rescaling_relative_error, arguments.termalization, particles,
      arguments.dt, L, LUT, arguments.verbose);
  } else {
    if (arguments.verbose) printf("No termalization, initial temperature will be the desired.\n");
    initial_velocities(particles, T);
  }

  update_forces(particles, L, LUT);

  // Evolve
  double time = 0;
  if (arguments.verbose) {
    printf("\n");
    printf("Beggining loop\n");
    printf("==============\n");
  }
  for (int frame = 0; frame < arguments.frames; frame++) {
    time = frame * arguments.frames_step * arguments.dt;
    if (arguments.verbose) printf(
      "\rAmount of frames covered: %d / %d (%.2lf%%)",
      frame, arguments.frames, 100 * (float)frame / arguments.frames);
    save_lammpstrj(arguments.lamppstrj_filename, particles.x, particles.v, N, L, frame);
    write_log(frame, time, arguments.log_filename, rho, L, LUT, particles);
    for (int i = 0; i < arguments.frames_step; i++) timestep(particles, arguments.dt, L, LUT);
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
