#ifndef lor_h
#define lor_h

#include "vector_ops.h"

typedef struct _lor {
  vec3* center;
  vec3* dir;
  double long_uncert;
  double transverse_uncert;
} lor;

typedef struct _rend render;

struct _rend {
  double* volume;     // the voxel volume
  int dimensions[3];  // number of voxels in each axis
  vec3* conversion;   // multiply by to convert from length to voxel space
  vec3* least_corner; // minimum position values of render volume
  vec3* max_corner;   // maximum position values of render volume
  void (*combiner)(render*, lor*); // adds a lor to the render
};

#endif