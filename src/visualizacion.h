#ifndef VISUALIZACION_H
#define VISUALIZACION_H

int save_lammpstrj(char *filename, float* x, float* v, int N, float L, int frame);
int load_frame(void *fp, float* x, float* v, int N, float *L);
int load_lammpstrj(char *filename, float* x, float* v, int N, float *L, int frame);

#endif
