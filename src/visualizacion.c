#include "visualizacion.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


int save_lammpstrj(char *filename, double* x, double* v, int N, double L, int frame){
  /* Save a frame on LAMMPSTRJ format.

  If it is the first frame, it initializes the file.
  */
  FILE *fp;

  if (frame) fp = fopen(filename, "a");
  else fp = fopen(filename, "w");

  // Header.
	fprintf(fp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n", frame, N);
	for(int l = 0; l < 3; l++) fprintf(fp, "0 %lf\n", L);

  // Atoms positions.
  fprintf(fp, "ITEM: ATOMS id x y z vx vy vz \n");
	for(int i = 0; i < N; i++) fprintf(
    fp, "%d %lf %lf %lf %lf %lf %lf\n",
    i,  // Particle number
    x[3 * i], x[3 * i + 1], x[3 * i + 2],  // Position for each particle.
    v[3 * i], v[3 * i + 1], v[3 * i + 2]  // Velocity for each particle.
  );

  fclose(fp);
  return 0;
}


int load_frame(void *fp, double* x, double* v, int N, double *L){
  /* Load LAMPPSTRJ frame into memory.

  Does so by reading through LAMMPSTRJ file stream. This way, if the
  stream has already been read, a new frame will be parsed. Given that
  the LAMMPSTRJ format is the following:

  ITEM: TIMESTEP
  {FRAME}
  ITEM: NUMBER OF ATOMS
  {N}
  ITEM: BOX BOUNDS pp pp pp
  0 {SIZE}
  0 {SIZE}
  0 {SIZE}
  ITEM: ATOMS id x y z vx vy vz
  {ID} {X} {Y} {Z} {VX} {VY} {VZ}
  {ID} {X} {Y} {Z} {VX} {VY} {VZ}
  ...

  Then one will read through the file and load this data accordingly.

  Returns
  -------
  * If -2, the file corresponds to a different amount of particles.
  * If -1, the file has no more content.
  * Else, the number of the frame loaded.
  */
  char buffer[255], *eof;
  int id, frame, N_file;

  // Get timestep (frame):
  eof = fgets(buffer, 255, fp);  // Gets first line (ITEM: TIMESTEP)
  if (eof == NULL) return -1;  // The file has no more information.
  id = fscanf(fp, "%d\n", &frame);  // Cargo el frame

  // Assert number of particles (N)
  fgets(buffer, 255, fp);  // Gets third line (ITEM: NUMBER OF ATOMS)
  id = fscanf(fp, "%d\n", &N_file); // Cargo el numero de particulas
  if (N_file != N) return -2;  // Has different amount of particles than stated

  // Get box size.
  for(int l = 0; l < 2; l++) eof = fgets(buffer, 255, fp);
  id = fscanf(fp, "0 %lf\n", L);  // Load box size
  for(int l = 0; l < 2; l++) eof = fgets(buffer, 255, fp);

  for(int i = 0; i < N; i++){  // Load positions and velocities of each particle
    id = fscanf(
      fp, "%d %lf %lf %lf %lf %lf %lf\n",
      &id,  // Id (unnecesary, as it is corresponded with position)
      x + 3 * i, x + 3 * i + 1, x + 3 * i + 2,  // Position (x, y ,z)
      v + 3 * i, v + 3 * i + 1, v + 3 * i + 2);  // Velocity (vx, vy ,vz)
  }
  return frame;
}


int load_lammpstrj(char *filename, double* x, double* v, int N, double *L, int frame){
  /* Load specific frame from LAMMPSTRJ filepath into x and v pointers.*/
  FILE *fp = fopen(filename, "r");
  int frame_file = load_frame(fp, x, v, N, L);
  while(frame_file < frame && frame_file >= 0){  // Load frames until frame is right
    frame_file = load_frame(fp, x, v, N, L);
    if(frame_file == -2) printf("Number of particles incompatible with file.\n");
  }
  fclose(fp);
  if(frame_file == -1) printf("Frame %d not in %s.\n", frame, filename);
  return frame_file;
}
