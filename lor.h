#ifndef lor_h
#define lor_h

#include "vector_ops.h"

typedef struct _lor {
  vec3* center;
  vec3* dir;
  double long_uncert;
  double transverse_uncert;
} lor;

#endif