#ifndef ray_trace_h
#define ray_trace_h

#include "vector_ops.h"

#define REC_PRISM 1
#define SPHERE 2
#define CYLINDER 3

typedef struct _shape_holder shape;
typedef struct _ray_trace ray;
typedef struct _traversal traversal;
typedef struct _multi_shape geometry;
typedef unsigned int uint;

struct _shape_holder {
  int type; // a holder for what type of 
  float pos[3]; // 3d location of the shape
  int axis; // direction that the primary axis of the shape is along (1 = x, 2 = y, 3 = z)
  float dim[3]; // dimentions needed to define the shape. dist across x,y,z, r, or r, h
  float atten; // the linear attenuation coefficent of the shape
};

struct _ray_trace {
  vec3d pos;
  vec3d dir;
};

struct _traversal {
  vec3d intersection;
  double t;
};

struct _multi_shape {
  shape** geo;
  uint size;
};

// builds a ray structure
ray* ray_build(vec3d pos, vec3d dir);

// frees a ray structure
void* ray_free(ray* src);

// frees a traversal structure
void traversal_free(traversal* src);

// makes a shape
shape* shape_build(int type, float* pos, float* dim, int axis, float attenuation);

// frees a shape (a free with return null)
void* shape_free(void* a);

// forms a geometry structure
geometry* geometry_build(shape** geos, uint size);

void* geometry_free(geometry* a);


/*
 * takes a ray and finds how far it spends inside of the given geometry. To do
 * this it steps through each object in the geometry. It checks to see if any
 * of the geometry has the path starting within the object. If it does, it marks
 * that geometry as crossed and starts again with the same requirements. If
 * the ray does not start in any of the geometry then the distance to all exits
 * is compared. The shortest distance from the ray start to the ray end less the
 * ray travel inside of the object is used as the next distance, and the
 * procedure continues from the exitpoint.
 */
double propagate(ray* path_src, geometry* all);

// the full suite of tests for ray tracing
int geometry_full_tests();

// prints a given geometry
void print_geometry(FILE* output, geometry* a);

/*
 * Takes a ray and prism and finds how far the ray traveled (in the direction of
 * travel) through the prism. Returns the exit location and distance travled.
 * If the traversal started outside the box (and therefore there are two crossings
 * in the direction of travel) the total distance inside the box is returned,
 * with the farther exit location and full_crossing is set to true.
 */
traversal* exit_rectangular_prism(ray* path, shape* prism, int* full_crossing);

/* 
 * Finds the distance to exiting a cylinder from the current ray. If the ray
 * starts outside and goes through it returns the distance traveled inside and
 * sets the full crossing flag to true. Works by checking for intersections with
 * the two flat ends, then looking for ones along the cylinder. The cylinder
 * check reuses the code for the sphere, just projected onto the xy plane.
 */
traversal* exit_cyl(ray* path, shape* cyl, int* full_crossing);

/*
 * finds the distance to exiting the sphere for the given ray. If the ray
 * starts outside the sphere and travels through, the distance is given as the
 * travel inside of the sphere and full_crossing is set to true. Works using
 * a line crossing sphere formula on wikipedia
 * 
 * Assumes that the ray is a unit vector
 */
traversal* exit_sphere(ray* path, shape* sphere, int* full_crossing);


#endif