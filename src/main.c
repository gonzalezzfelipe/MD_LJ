#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "aphrodite.h"  // Visualization
#include "arguments.h"  // Argument parsing
#include "chronos.h"    // Temporal evolution
#include "force.h"      // Forces interaction
#include "helpers.h"    // Helpers for improved experience
#include "init.h"       // Initialization of everything
#include "metrics.h"    // Important metrics calculation
#include "objects.h"    // Definition of objects used in program

// Default options
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
  strcpy(arguments.lammpstrj_filename, DEF_LAMMPSTRJ_FILENAME);
  strcpy(arguments.log_filename, DEF_LOG_FILENAME);
  arguments.seed = time(NULL);
  arguments.r_c = DEF_R_C;
  arguments.table_length = DEF_TABLE_LENGTH;
  arguments.rescaling_relative_error = DEF_RESCALING_RELATIVE_ERROR;
  arguments.initial_temperature = DEF_INITIAL_TEMPERATURE;
  arguments.verbose = 0;
  arguments.exact = 0;

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
    printf("Name of LAMPPSTRJ file: %s\n", arguments.lammpstrj_filename);
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

  if (!(arguments.exact)) fill_lut(&LUT, arguments.r_c, arguments.table_length);

  particles.x = (double*)malloc(3 * N * sizeof(double));
  particles.v = (double*)malloc(3 * N * sizeof(double));
  particles.f = (double*)malloc(3 * N * sizeof(double));
  particles.N = N;

  initial_positions(L, particles);
  update_forces(particles, L, LUT, arguments.exact);

  // Termalize
  if (arguments.termalization) {
    initial_velocities(particles, arguments.initial_temperature);
    rescaling(
      T, arguments.rescaling_relative_error, arguments.termalization, particles,
      arguments.dt, L, LUT, arguments.verbose, arguments.exact);
  } else {
    if (arguments.verbose) printf("No termalization, initial temperature will be the desired.\n");
    initial_velocities(particles, T);
  }

  // Evolve
  double time = 0;
  if (arguments.verbose) printf("\nBeggining loop\n==============\n");
  for (int frame = 0; frame < arguments.frames; frame++) {
    time = frame * arguments.frames_step * arguments.dt;
    if (arguments.verbose) progress(frame, arguments.frames);
    save_lammpstrj(arguments.lammpstrj_filename, particles.x, particles.v, N, L, frame);
    write_log(
      frame, time, arguments.log_filename, rho, L, LUT, particles, arguments.exact);
    for (int i = 0; i < arguments.frames_step; i++)
      timestep(particles, arguments.dt, L, LUT, arguments.exact);
  }
  progress(arguments.frames, arguments.frames);

  free(particles.x);
  free(particles.v);
  free(particles.f);
  free(LUT.r);
  free(LUT.r2);
  free(LUT.v);
  free(LUT.f);
  return 0;
}
