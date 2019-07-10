#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define _EMPTY '_'
#define _FULL '#'
#define _TOTAL 70


unsigned long get_epoch_in_microseconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long ret = tv.tv_usec;
  ret += (tv.tv_sec * 1000000);
  return ret;
}


int progress(int done, int total, unsigned long start) {
  /* Print progressbar to console. */
  char _progress[_TOTAL];

  float lapsed = (float)(get_epoch_in_microseconds() - start);
  int eta = (int)(lapsed * (total / 1.0 / done - 1) / 1000000);
  if (eta < 0) eta = 0;

  int covered = (int)(_TOTAL * done / 1.0 / total);
  for (int i = 0; i < _TOTAL; i++) {
    if (i < covered) _progress[i] = _FULL;
    else _progress[i] = _EMPTY;
  }
  printf(
    "\rProgress: |%s| (%.2f %%) ETA: %02d:%02d:%02d ",
    _progress, 100 * done / 1.0 / total,
    eta / 3600,  // Hours
    (eta % 3600) / 60,  // Minutes
    eta % 60  // Seconds
  );
  fflush(stdout);
  if (done == total) printf("\n");
  return 0;
}


#endif
