#ifndef OBJETOS_H
#define OBJETOS_H

typedef struct Particles {
  double* x;
  double* v;
  double* f;
  int N;
} Particles;

typedef struct LookUpTable {
  double* r;
  double* r2;
  double* v;
  double* f;
  double r_c;
  int length;
} LookUpTable;

#endif
