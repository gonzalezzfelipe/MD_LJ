#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>

#define _EMPTY '_'
#define _FULL '#'
#define _TOTAL 50


int progress(int done, int total) {
  /* Print progressbar to console. */
  char _progress[_TOTAL];

  int covered = (int)(_TOTAL * done / 1.0 / total);
  for (int i = 0; i < _TOTAL; i++) {
    if (i < covered) _progress[i] = _FULL;
    else _progress[i] = _EMPTY;
  }
  printf("\rProgress: %s (%.2f %%)", _progress, 100 * done / 1.0 / total);
  fflush(stdout);
  if (done == total) printf("\n");
  return 0;
}


#endif
