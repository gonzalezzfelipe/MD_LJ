/* Module for parsing command line arguments. */
#ifndef ARGUMENTOS_H
#define ARGUMENTOS_H

#define AMOUNT_OF_POSITIONAL_ARGS 3

#include <argp.h>
#include <stdlib.h>
#include <string.h>

// Version
const char *argp_program_version = "md 0.1.0";
const char *argp_program_bug_address = "<gonzalezz_felipe@hotmail.com>";

// Documentation of the program.
static char doc[] =
  "Program to simulate molecular dynamics on a periodic cubic box.";

// Required, positional arguments
static char args_doc[] = "NUMBER_OF_PARTICLES DENSITY TEMPERATURE";

// Keyword arguments.
static struct argp_option options[] = {
  {"dt", 'h', "DT", 0, "Temporal step." },
  {"frames", 'f', "FRAMES", 0, "Amount of frames to save on output files."},
  {"frames_step", 's', "FRAMES_STEP", 0, "Amount of frames to leave between saved frames."},
  {"termalization", 't', "TERMALIZATION", 0, "Amount of timesteps for termalization and rescaling."},
  {"lammpstrj_filename", 'o', "LAMMPSTRJ_FILENAME", 0, "Filename where to write LAMMPS trajectories."},
  {"log_filename", 'l', "LOG_FILENAME", 0, "Filename where to log relevant metrics."},
  {"seed", 'r', "SEED", 0, "Random seed."},
  {"r_c", 'c', "CUT_RADIUS", 0, "Radius for which interaction becomes 0."},
  {"table_length", 'q', "LUT_LENGTH", 0, "Length of LUT table."},
  {"rescaling_relative_error", 'w', "RESCALING_RELATIVE_ERROR", 0, "Relative error for temperature rescaling."},
  {"initial_temperature", 'i', "INITIAL_TEMPERATURE", 0, "Initial temperature before rescaling."},
  {"verbose", 'v', 0, 0, "Whether to be verbose on the output."},
  {"exact", 'e', 0, 0, "If defined, table will not be used and the forces and potentials will be calculated."},
  { 0 }
};

// Struct where arguments will be stored.
struct arguments {
  // First the positional args
  char *args[AMOUNT_OF_POSITIONAL_ARGS];

  // Options
  double dt;
  double r_c;
  int table_length;
  double rescaling_relative_error;
  double initial_temperature;
  int frames;
  int frames_step;
  int termalization;
  char lammpstrj_filename[255];
  char log_filename[255];
  int seed;

  // Flags
  int verbose;
  int exact;
};


// Parser that should handle a single option each time.
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Parse a argument.

  For any program a similar structure can be defined. Each key has to have
  a case in which the argument is correctly parsed.
  */
  struct arguments *arguments = state->input;

  switch (key) {
    case 'h':
      sscanf(arg, "%lf", &arguments->dt);
      break;
    case 'f':
      sscanf(arg, "%d", &arguments->frames);
      break;
    case 's':
      sscanf(arg, "%d", &arguments->frames_step);
      break;
    case 't':
      sscanf(arg, "%d", &arguments->termalization);
      break;
    case 'o':
      strcpy(arguments->lammpstrj_filename, arg);
      break;
    case 'l':
      strcpy(arguments->log_filename, arg);
      break;
    case 'r':
      sscanf(arg, "%d", &arguments->seed);
      break;
    case 'c':
      sscanf(arg, "%lf", &arguments->r_c);
      break;
    case 'q':
      sscanf(arg, "%d", &arguments->table_length);
      break;
    case 'w':
      sscanf(arg, "%lf", &arguments->rescaling_relative_error);
      break;
    case 'i':
      sscanf(arg, "%lf", &arguments->initial_temperature);
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'e':
      arguments->exact = 1;
      break;

    // This below is general to any program.
    case ARGP_KEY_ARG:
      if (state->arg_num >= AMOUNT_OF_POSITIONAL_ARGS) argp_usage(state);
      arguments->args[state->arg_num] = arg;
      break;

    case ARGP_KEY_END:
      if (state->arg_num < AMOUNT_OF_POSITIONAL_ARGS) argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}


// The actual argument parser
static struct argp argp = {options, parse_opt, args_doc, doc};


#endif
