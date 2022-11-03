#ifndef lor_h
#define lor_h

#include "vector_ops.h"

typedef unsigned int uint;

typedef struct _lor {
  vec3d* center;
  vec3d* dir;
  double long_uncert;
  double transverse_uncert;
} lor;

typedef struct _rend render;

struct _rend {
  double* volume;     // the voxel volume
  int dimensions[3];  // number of voxels in each axis
  vec3d* conversion;   // multiply by to convert from length to voxel space
  vec3d* least_corner; // minimum position values of render volume
  vec3d* max_corner;   // maximum position values of render volume
  double cutoff;      // the cutoff beyond which no values are processed for a LOR
  void (*combiner)(render*, double, int*, double); // adds a value of the lor to the render
};

render* master_copy = NULL;

#endif