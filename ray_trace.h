#ifndef ray_trace_h
#define ray_trace_h

#include "vector_ops.h"

#define REC_PRISM 1
#define SPHERE 2
#define CYLINDER 3

typedef struct _shape_holder shape;
typedef struct _ray_trace ray;
typedef struct _traversal traversal;

struct _shape_holder {
  int type; // a holder for what type of 
  float position[3]; // 3d location of the shape
  int axis; // direction that the primary axis of the shape is along (1 = x, 2 = y, 3 = z)
  float dimentions[3]; // dimentions needed to define the shape. dist across x,y,z
  float attenuation; // the linear attenuation coefficent of the shape
  vec3d* (*intersection)(shape*, ray*); // function that provides when a ray
  // crosses the surface of a shape
};

struct _ray_trace {
  vec3d* pos;
  vec3d* dir;
};

struct _traversal {
  vec3d* intersection;
  double t;
};

#endif